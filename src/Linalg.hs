{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}

import Data.Matrix
import Data.Bifunctor
import qualified Data.Vector as V
import Debug.Trace as DT
import Data.CReal

subtr :: Num a => Matrix a -> Matrix a -> Matrix a
subtr = elementwise (-)
 
add :: Num a => Matrix a -> Matrix a -> Matrix a
add = elementwise (+)

norm2 :: Floating a => Matrix a -> a
norm2 v = sqrt $ foldr (\x r -> x * x + r) 0 (toList v) 

getDiagonal :: Num a => Matrix a -> Matrix a
getDiagonal m = diagonalList (nrows m) 0 $ V.toList $ getDiag m

gershgorinCircles :: Floating a => Matrix a -> [(a, a)]
gershgorinCircles m = zip cs rs
    where
        cs   = V.toList $ getDiag m
        diag = getDiagonal m
        rs   = map (foldr (\x r -> abs x + r) 0) (toLists $ m `subtr` diag)

outUnitCircle :: (Floating a, Ord a) => Matrix a -> Bool
outUnitCircle m = foldr (\(c, r) acc -> (abs c + r >= 1) || acc) False (gershgorinCircles m)

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
    | norm2 (x `subtr` x') <= eps    = Right x'
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
    | norm2 ((l `multStd` x) `subtr` b') <= eps = Right x
    | otherwise                                 = doGaussZeidel l negu b x' eps outCircle 0
    where
        b' = (negu `multStd` x) `add` b
        x' = gaussPartial l b'

givensRotation :: Floating a => [[a]] -> Int -> Int -> a -> a -> [[a]]
givensRotation m i j c s = zipWith f [1..] m
    where
        ui = m !! (i - 1)
        uj = m !! (j - 1)
        uinew = zipWith (\xi xj -> c * xi + s * xj) ui uj
        ujnew = zipWith (\xi xj -> (- s) * xi + c * xj) ui uj
        f pos row
            | pos == i  = uinew
            | pos == j  = ujnew
            | otherwise = row

trans :: [[a]] -> [[a]]
trans = toLists . transpose . fromLists

idMatrix :: Floating a => Int -> [[a]]
idMatrix sz = toLists $ identity sz

firstNonzero :: (Floating a, Ord a) => [a] -> Int -> Int
firstNonzero v k = foldl (\acc (j, x) -> if j >= k && x /= 0 && acc == 0 
                                       then j 
                                       else acc) 0 (zip [1..] v) 

qrDecompGivens :: (Floating a, Ord a) => [[a]] -> (Matrix a, Matrix a)
qrDecompGivens m = appToPair fromLists $ first trans qr
    where
        appToPair f (x, y) = (f x, f y)
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
                    | j <= k    = ((q', r'), xi)
                    | otherwise = ((givensRotation q' j i c s, givensRotation r' j i c s), (- s) * xj + c * xi)
                        where
                            n = sqrt $ xi * xi + xj * xj
                            c = xi / n
                            s = (- xj) / n

multHouseholder :: Floating a => [[a]] -> [[a]] -> [[a]]
multHouseholder m' v' = toLists $ m `subtr` prod
    where
        m   = fromLists m'
        v   = fromLists v'
        vt  = transpose v
        prod = scaleMatrix 2 v `multStd` (vt `multStd` m)

qrDecompHouseholder :: (Floating a, Ord a) => [[a]] -> (Matrix a, Matrix a)
qrDecompHouseholder m = appToPair fromLists $ first trans qr
    where
        appToPair f (x, y) = (f x, f y)
        idm                = idMatrix $ length m
        qr                 = foldl handler (idm, m) [1..(length m - 1)]
        handler (q, r) k 
            | i == 0 || u == e1 = (q, r)
            | otherwise = (multHouseholder q v, multHouseholder r v)
            where
                col = trans r !! (k - 1)
                i   = firstNonzero col k
                v'  = fromLists $ zipWith (\j x -> if j < i then [0] else [x]) [1..length col] col
                u   = scaleMatrix (1 / norm2 v') v'
                e1  = mapPos (\(j, _) _-> if j == i then 1 else 0) v'
                v   = toLists $ scaleMatrix (1 / norm2 (u `subtr` e1)) (u `subtr` e1)