{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE InstanceSigs #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE UndecidableInstances #-}
 
import Data.Matrix
import qualified Data.Vector
 
subtr :: Num a => Matrix a -> Matrix a -> Matrix a
subtr = elementwise (-)
 
add :: Num a => Matrix a -> Matrix a -> Matrix a
add = elementwise (+)
 
norm2 :: Floating a => Matrix a -> a
norm2 v = sqrt $ foldr (\x r -> x * x + r) 0 (toList v) 

getDiagonal :: (Floating a, Ord a) => Matrix a -> Matrix a
getDiagonal m = diagonalList (nrows m) 0 $ Data.Vector.toList $ getDiag m

gershgorinCircles :: (Floating a, Ord a) => Matrix a -> [(a, a)]
gershgorinCircles m = zip cs rs
    where
        cs   = Data.Vector.toList $ getDiag m
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
 
lowerTriangle :: (Floating a, Ord a) => Matrix a -> Matrix a
lowerTriangle m = mapPos (\(i, j) e -> if (size - i + 1) + j <= size + 1 then e else 0) m
    where
        size = nrows m 

upperTriangle :: (Floating a, Ord a) => Matrix a -> Matrix a
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