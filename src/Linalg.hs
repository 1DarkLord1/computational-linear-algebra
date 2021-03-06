module Linalg where

{-# LANGUAGE RankNTypes,
             FlexibleContexts,
             FlexibleInstances,
             MultiParamTypeClasses,
             UndecidableInstances,
             AllowAmbiguousTypes,
             DataKinds,
             ScopedTypeVariables,
#-}


import Data.Matrix
import Data.Bifunctor
import qualified Data.Vector as V
import qualified Data.List.Zipper as ZP
import Data.Complex
import qualified Data.List as DL
import qualified Data.Set as DS
import qualified Data.Map as DM
import GHC.Float
import Data.Maybe
import Data.Euclidean

subtr :: Num a => Matrix a -> Matrix a -> Matrix a
subtr = elementwise (-)
 
add :: Num a => Matrix a -> Matrix a -> Matrix a
add = elementwise (+)

norm2 :: Floating a => Matrix a -> a
norm2 v = sqrt $ foldr (\x r -> x * x + r) 0 (toList v) 

norm2C :: RealFloat a => Matrix (Complex a) -> a
norm2C v = sqrt $ foldr (\x r -> realPart (x * conjugate x) + r) 0 (toList v) 

getDiagonal :: Num a => Matrix a -> Matrix a
getDiagonal m = diagonalList (nrows m) 0 $ V.toList $ getDiag m

gershgorinCircles :: Num a => Matrix a -> [(a, a)]
gershgorinCircles m = zip cs rs
    where
        cs   = V.toList $ getDiag m
        diag = getDiagonal m
        rs   = map (foldr (\x r -> abs x + r) 0) (toLists $ m `subtr` diag)

outUnitCircle :: (Num a, Ord a) => Matrix a -> Bool
outUnitCircle m = foldr (\(c, r) acc -> (abs c + r >= 1) || acc) False (gershgorinCircles m)


{- Simple iteration method of solving linear equations systems. Complexity O(n^2) per iteration. -}
simpleIteration :: (Floating a, Ord a) => [[a]] -> [[a]] -> a -> Either Int (Matrix a)
simpleIteration m b = simpleIteration' (fromLists m) (fromLists b)
  
simpleIteration' :: (Floating a, Ord a) => Matrix a -> Matrix a -> a -> Either Int (Matrix a)
simpleIteration' m b eps = doIterations m' b (zero size 1) eps (outUnitCircle m) 0 
    where 
          size = nrows m
          m'   = identity size `subtr` m
 
doIterations m b x eps outCircle cnt
    | outCircle 
        && cnt == 20                 = Left 0
    | outCircle
        && norm2 x' >= norm2 x + 1   = doIterations m b x' eps outCircle (cnt + 1) 
    | norm2 (x `subtr` x') < eps     = Right x'
    | otherwise                      = doIterations m b x' eps outCircle 0
    where x'                         = (m `multStd` x) `add` b
 
lowerTriangle :: Num a => Matrix a -> Matrix a
lowerTriangle m = mapPos (\(i, j) e -> if (size - i + 1) + j <= size + 1 then e else 0) m
    where
        size = nrows m 

upperTriangle :: Num a => Matrix a -> Matrix a
upperTriangle m = transpose $ m' `subtr` diag 
    where
        m'   = lowerTriangle $ transpose m  
        diag = getDiagonal m

gaussPartial :: (Floating a, Ord a) => Matrix a -> Matrix a -> Matrix a
gaussPartial m b = fromLists $ map (: []) (foldl gaussStep [] mlist) 
    where
        m'                = b <|> m
        mlist             = toLists m'
        gaussStep ans row = ans ++ [xi]
            where
                coef = foldr (\e acc -> if e /= 0 && acc == 0 then e else acc) 0 row
                bi   = head row
                row' = tail row
                sm   = sum $ zipWith (*) ans row'
                xi   = (bi - sm) / coef  

{- Gauss-Zeidel method of solving linear equations systems. Complexity O(n^2) per iteration. -}
gaussZeidel :: (Floating a, Ord a) => [[a]] -> [[a]] -> a -> Either Int (Matrix a)
gaussZeidel m b = gaussZeidel' (fromLists m) (fromLists b)

gaussZeidel' :: (Floating a, Ord a) => Matrix a -> Matrix a -> a -> Either Int (Matrix a)
gaussZeidel' m b eps = doGaussZeidel l negu b (zero size 1) eps (outUnitCircle m) 0 
    where
        size = nrows m
        l    = lowerTriangle m
        negu = zero size size `subtr` upperTriangle m

doGaussZeidel l negu b x eps outCircle cnt
    | outCircle 
        && cnt == 20                            = Left 0
    | outCircle
        && norm2 x' >= norm2 x + 1              = doGaussZeidel l negu b x' eps outCircle (cnt + 1)
    | norm2 ((l `multStd` x) `subtr` b') < eps  = Right x
    | otherwise                                 = doGaussZeidel l negu b x' eps outCircle 0
    where
        b' = (negu `multStd` x) `add` b
        x' = gaussPartial l b'

{- Givens rotation of matrix. Complexity O(n). -}
givensRotation :: Num a => [[a]] -> Int -> Int -> a -> a -> [[a]]
givensRotation m i j c s = zipWith subst [1..] m
    where
        ui    = m !! (i - 1)
        uj    = m !! (j - 1)
        uinew = zipWith (\xi xj -> c * xi + s * xj) ui uj
        ujnew = zipWith (\xi xj -> (- s) * xi + c * xj) ui uj
        subst pos row
            | pos == i  = uinew
            | pos == j  = ujnew
            | otherwise = row

trans :: [[a]] -> [[a]]
trans = toLists . transpose . fromLists

idMatrix :: Num a => Int -> [[a]]
idMatrix sz = toLists $ identity sz

firstNonzero :: (Num a, Ord a) => [a] -> Int -> Int
firstNonzero v k = foldl (\acc (j, x) -> if j >= k && x /= 0 && acc == 0 
                                       then j 
                                       else acc) 0 (zip [1..] v) 

appToPair f (x, y) = (f x, f y)

{- Givens rotation method of QR decomposition. Complexity O(n^3). -}
qrDecompGivens :: (Floating a, Ord a) => [[a]] -> (Matrix a, Matrix a)
qrDecompGivens m = appToPair fromLists $ first trans qr
    where
        idm                = idMatrix $ length m 
        qr                 = foldl handler (idm, m) [1..(length m - 1)]
        handler (q, r) k 
            | i == 0    = (q, r)
            | otherwise = (givensRotation q'' k i 0 1, givensRotation r'' k i 0 1)
            where
                col        = trans r !! (k - 1)
                i          = firstNonzero col k
                (q'', r'') = fst $ foldl handler' ((q, r), col !! (i - 1)) (zip [1..] col)
                handler' ((q', r'), xi) (j, xj)
                    | j <= i    = ((q', r'), xi)
                    | otherwise = ((givensRotation q' j i c s, givensRotation r' j i c s), (- s) * xj + c * xi)
                        where
                            n = sqrt $ xi * xi + xj * xj
                            c = xi / n
                            s = (- xj) / n

{- Householder matrix multiplication. Complexity O(n^2). -}
multHouseholder :: Num a => [[a]] -> [[a]] -> [[a]]
multHouseholder m' v' = toLists $ m `subtr` prod
    where
        m    = fromLists m'
        v    = fromLists v'
        vt   = transpose v
        prod = scaleMatrix 2 v `multStd` (vt `multStd` m)

{- Householder reflection method of QR decomposition. Complexity O(n^3). -}
qrDecompHouseholder :: (Floating a, Ord a) => [[a]] -> (Matrix a, Matrix a)
qrDecompHouseholder m = appToPair fromLists $ first trans qr
    where
        idm                = idMatrix $ length m
        qr                 = foldl handler (idm, m) [1..(length m - 1)]
        handler (q, r) k 
            | norm2 v' == 0 || u == e1 = (q, r)
            | otherwise                = (multHouseholder q v, multHouseholder r v)
            where
                col = trans r !! (k - 1)
                v'  = fromLists $ zipWith (\j x -> if j < k then [0] else [x]) [1..length col] col
                u   = scaleMatrix (1 / norm2 v') v'
                e1  = mapPos (\(j, _) _-> if j == k then 1 else 0) v'
                v   = toLists $ scaleMatrix (1 / norm2 (u `subtr` e1)) (u `subtr` e1)
 
{- Simple iteration method of calculation the maximum modulo eigenvalue. Complexity O(n^2) per iteration. -}
simpleIterationMaxEV :: RealFloat a => [[Complex a]] -> [[Complex a]] -> a -> Int
                        -> Either Int (Complex a, Matrix (Complex a))
simpleIterationMaxEV m v = doItersMaxEV (fromLists m) (fromLists v)

doItersMaxEV m x eps cnt
    | cnt == 0          = Left 0
    | norm2C diff < eps = Right (ev, x)
    | otherwise         = doItersMaxEV m x' eps (cnt - 1) 
    where
        norm = norm2C (m `multStd` x) :+ 0
        x'   = scaleMatrix (1 / norm) (m `multStd` x)
        ev   = head $ head $ toLists $ transpose x `multStd` (m `multStd` x)
        diff = (m `multStd` x) `subtr` scaleMatrix ev x

{- QR-algorithm of calculation the matrix spectrum. Complexity O(n^3) per iteration. -}
qrEV :: (Floating a, Ord a) => [[a]] -> a -> ([a], Matrix a)
qrEV m eps = doItersQrEV m eps (identity $ length m)

doItersQrEV m eps qk
    | lessEps   = (evs, qk)
    | otherwise = doItersQrEV m' eps (qk `multStd` q)
    where
        (q, r)  = qrDecompGivens m
        m'      = toLists $ r `multStd` q
        circles = gershgorinCircles (fromLists m)
        lessEps = foldr (\(_, rd) acc -> (rd < eps) && acc) True circles
        evs     = V.toList $ getDiag $ fromLists m

multHouseholderRight :: Num a => [[a]] -> [[a]] -> [[a]]
multHouseholderRight m' v' = toLists $ m `subtr` prod
    where
        m    = fromLists m'
        v    = fromLists v'
        vt   = transpose v
        prod = scaleMatrix 2 ((m `multStd` v) `multStd` vt)

{- Matrix transformation to tridiagonal form. Complexity O(n^3). -}
getTridiagonal :: (Floating a, Ord a) => [[a]] -> (Matrix a, Matrix a)
getTridiagonal m = appToPair fromLists $ second trans tridiag
    where
        idm                = idMatrix $ length m
        tridiag            = foldl handler (m, idm) [1..(length m - 1)]
        handler (a, q) k 
            | norm2 v' == 0 || u == e1 = (a, q)
            | otherwise                = (newa, newq)
            where
                col  = trans a !! (k - 1)
                v'   = fromLists $ zipWith (\j x -> if j <= k then [0] else [x]) [1..length col] col
                u    = scaleMatrix (1 / norm2 v') v'
                e1   = mapPos (\(j, _) _-> if j == k + 1 then 1 else 0) v'
                v    = toLists $ scaleMatrix (1 / norm2 (u `subtr` e1)) (u `subtr` e1)
                newa = multHouseholderRight (multHouseholder a v) v
                newq = multHouseholder q v

cursorp :: ZP.Zipper a -> Int
cursorp (ZP.Zip l _) = length l

rightn :: Int -> ZP.Zipper a -> ZP.Zipper a
rightn 0 z = z
rightn n z = rightn (n - 1) (ZP.right z) 

givensRotationZ :: Num a => [ZP.Zipper a] -> Int -> Int -> a -> a -> [ZP.Zipper a]
givensRotationZ m i j c s = zipWith subst [1..] m
    where
        ui      = m !! (i - 1)
        uj      = m !! (j - 1)
        curspUi = cursorp ui
        curspUj = cursorp uj   
        uinew   = rightn curspUi $ ZP.fromList $ 
                   zipWith (\xi xj -> c * xi + s * xj) (ZP.toList ui) (ZP.toList uj)
        ujnew   = rightn curspUj $ ZP.fromList $ 
                   zipWith (\xi xj -> (- s) * xi + c * xj) (ZP.toList ui) (ZP.toList uj)
        subst pos row
            | pos == i  = uinew
            | pos == j  = ujnew
            | otherwise = row

{- QR decomposition for the tridiagonal matrices. Complexity O(n^2). -}
qrDecompTridiagonal :: (Floating a, Ord a) => [[a]] -> ([(Int, Int, a, a)], Matrix a)
qrDecompTridiagonal m = second (fromLists . fmap ZP.toList) $ 
                        foldl handler ([], zipped_m) [1..(length m - 1)]
    where
        zipped_m = fmap ZP.fromList m
        handler (gs, r) k
            | n == 0    = (gs, r)
            | otherwise = ((k + 1, k, c, s) : gs, map ZP.right $ givensRotationZ r (k + 1) k c s)
            where
                col = map ZP.cursor r
                xi  = col !! (k - 1)
                xj  = col !! k
                n   = sqrt $ xi * xi + xj * xj
                c   = xi / n
                s   = (- xj) / n

multGivens :: Num a => [(Int, Int, a, a)] -> [[a]] -> [[a]]
multGivens gs m = foldr (\(i, j, c, s) acc -> givensRotation acc i j c s) m gs

{- 
QR-algorithm of calculation the matrix spectrum for the tridiagonal matrices. 
Complexity O(n^2) per iteration. 
-}
qrEVTridiagonal :: (Floating a, Ord a) => [[a]] -> a -> ([a], Matrix a)
qrEVTridiagonal m eps = second (transpose . fromLists) $ doItersQrEVsTridiagonal m eps (idMatrix $ length m)

doItersQrEVsTridiagonal m eps qk
    | lessEps   = (evs, qk)
    | otherwise = doItersQrEVsTridiagonal m' eps q
    where
        (gs, r) = qrDecompTridiagonal m
        m'      = trans $ multGivens gs (trans $ toLists r)
        q       = multGivens gs qk
        circles = gershgorinCircles (fromLists m)
        lessEps = foldr (\(_, rd) acc -> (rd < eps) && acc) True circles
        evs     = V.toList $ getDiag $ fromLists m

doItersQrEVsTridiagonalNTimes m cnt qk
    | cnt == 0  = (evs, qk)
    | otherwise = doItersQrEVsTridiagonalNTimes m' (cnt - 1) q
    where
        (gs, r) = qrDecompTridiagonal m
        m'      = trans $ multGivens gs (trans $ toLists r)
        q       = multGivens gs qk
        evs     = V.toList $ getDiag $ fromLists m

doItersQrMinEVTridiagonal m eps qk
    | radius < eps = (qk, m)
    | otherwise    = doItersQrMinEVTridiagonal m' eps q
    where
        (gs, r) = qrDecompTridiagonal m
        m'      = trans $ multGivens gs (trans $ toLists r)
        q       = multGivens gs qk
        lastRow = last m
        radius  = foldr (\x acc -> abs x + acc) 0 lastRow - abs (last lastRow)
 
swapMinor :: Num a => Matrix a -> Matrix a -> Matrix a
swapMinor minor m = m' `add` minor' 
    where
        minorSz = nrows minor
        mSz     = nrows m    
        m'     = mapPos (\(i, j) x -> if i <= minorSz && j <= minorSz
                                  then 0
                                  else x) m
        minor' = extendTo 0 mSz mSz minor

{- 
Wilkinson shifts QR-algorithm of calculation the matrix spectrum for the tridiagonal matrices. 
Complexity O(n^2) per iteration.
-}
qrEVShifts :: (Floating a, Ord a) => [[a]] -> a -> ([a], Matrix a)
qrEVShifts mt eps = (evs, transpose q)
    where
        sz     = length mt
        (q, d) = foldl handler (identity (length mt), fromLists mt) [0..sz - 2]
        evs    = V.toList $ getDiag d
        handler (qk, m) k 
               = (qt `multStd` qk, swapMinor r' m)
               where
                    m'       = submatrix 1 (sz - k) 1 (sz - k) m
                    sz'      = nrows m'
                    lowerSq  = submatrix (sz' - 1) sz' (sz' - 1) sz' m'
                    maxIters = 20
                    s        = head $ tail $ fst $ 
                               doItersQrEVsTridiagonalNTimes (toLists lowerSq) maxIters (idMatrix 2)
                    m''      = m' `subtr` scaleMatrix s (identity sz')
                    (qt, r)  = first fromLists $ doItersQrMinEVTridiagonal (toLists m'') eps (idMatrix sz)
                    r'       = fromLists r `add` scaleMatrix s (identity sz')

{- Test graphs on nonisomorphism. Complexity O(Time(qrEVShifts)). -}
isomorphic :: (Floating a, Ord a) => [[a]] -> [[a]] -> Bool
isomorphic g1 g2
    | length g1 /= length g2 = False
    | otherwise              = norm2 (g1Spec `subtr` g2Spec) <= 10 * eps
    where
        eps    = 0.00001
        g1'    = toLists $ fst $ getTridiagonal g1
        g2'    = toLists $ fst $ getTridiagonal g2
        g1Spec = fromLists $ map (: []) $ DL.sort $ fst $ qrEVShifts g1' eps
        g2Spec = fromLists $ map (: []) $ DL.sort $ fst $ qrEVShifts g2' eps

getAdjMap :: Num a => Int -> DM.Map Int (DM.Map Int a)
getAdjMap n = adjMap 
    where
        buildRow = DM.fromList $ map (\k -> (k, 0)) [0..n - 1]
        adjMap   = DM.fromList $ map (\k -> (k, buildRow)) [0..n - 1]

addEdge (x, y) = DM.update (Just . DM.update (Just . (+ 1)) y) x

buildAdjMatrix :: Num a => Int -> (Int -> DM.Map Int (DM.Map Int a) -> DM.Map Int (DM.Map Int a)) -> Matrix a
buildAdjMatrix n genEdges = m
    where
        adjMap      = getAdjMap n
        adjMap'     = foldr genEdges adjMap [0..n - 1]
        m           = matrix n n fill
        fill (i, j) = x
            where
                i' = i - 1
                j' = j - 1
                x  = fromJust $ DM.lookup (j - 1) $ fromJust $ DM.lookup (i - 1) adjMap'

buildGraph1 :: Int -> Matrix Double
buildGraph1 n = buildAdjMatrix (n * n) genEdges
    where
        genEdges x = foldr (\y acc -> genEdges' (x, y) . acc) id [0..n - 1]
        genEdges' (x, y)
            | x >= n    = id
            | otherwise = addEdge (v, u1) . addEdge (v, u2) . 
                          addEdge (v, u3) . addEdge (v, u4) .
                          addEdge (v, u5) . addEdge (v, u6) . 
                          addEdge (v, u7) . addEdge (v, u8)
            where
                x1 = (x + 2 * y) `mod` n
                x2 = (x - 2 * y + 3 * n) `mod` n
                x3 = (x + 2 * y + 1) `mod` n
                x4 = (x - 2 * y - 1 + 3 * n) `mod` n
                y1 = (y + 2 * x) `mod` n
                y2 = (y - 2 * x + 3 * n) `mod` n
                y3 = (y + 2 * x + 1) `mod` n
                y4 = (y - 2 * x - 1 + 3 * n) `mod` n
                u1 = x1 * n + y
                u2 = x2 * n + y
                u3 = x3 * n + y
                u4 = x4 * n + y
                u5 = x * n + y1
                u6 = x * n + y2
                u7 = x * n + y3
                u8 = x * n + y4
                v  = x * n + y
                
inv :: Int -> Int -> Int
inv x md = ((xinv `mod` md) + md) `mod` md
    where 
        xinv = snd $ gcdExt x md  

buildGraph2 :: Int -> Matrix Double
buildGraph2 p = buildAdjMatrix (p + 1) genEdges
    where
        genEdges x  
            | x == p    = addEdge (x, xInv) . addEdge (x, x) . addEdge (x, x)
            | otherwise = addEdge (x, xInv) . addEdge (x, x1) . addEdge (x, x2)
            where
                xInv
                    | x == 0    = p
                    | x == p    = 0
                    | otherwise = inv x p
                x1 = (x + 1) `mod` p
                x2 = (x - 1 + p) `mod` p

{- Calculation the optimal alpha for expander. Complexity O(Time(qrEVShifts)). -}
expanderAlpha :: Int -> Matrix Double -> Double -> Double
expanderAlpha n g d = max (abs ev1) (abs ev2) / d 
    where
        eps       = 0.00001
        evs       = fst $ qrEVShifts (toLists $ fst $ getTridiagonal $ toLists g) eps
        sortedEVs = reverse $ DL.sort evs
        ev1       = head $ tail sortedEVs
        ev2       = last sortedEVs

expanderAlpha1 :: Int -> Double 
expanderAlpha1 n = expanderAlpha n g 8
    where
        g = buildGraph1 n

expanderAlpha2 :: Int -> Double
expanderAlpha2 n = expanderAlpha (n + 1) g 3
    where
        g = buildGraph2 n 