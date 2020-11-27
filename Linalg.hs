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
 
gershgorinCircles :: (Floating a, Ord a) => Matrix a -> [(a, a)]
gershgorinCircles m = zip cs rs
    where
        cs   = Data.Vector.toList $ getDiag m
        diag = diagonalList (nrows m) 0 cs
        rs   = map (foldr (\x r -> abs x + r) 0) (toLists $ m `subtr` diag)

simpleIteration :: (Floating a, Ord a) => [[a]] -> [[a]] -> a -> Either Int (Matrix a)
simpleIteration m b = simpleIteration' (fromLists m) (fromLists b)
 
simpleIteration' :: (Floating a, Ord a) => Matrix a -> Matrix a -> a -> Either Int (Matrix a)
simpleIteration' m b eps = doIterations m' b (fromList size 1 zero) eps outUnitCircle 0 
    where 
          m' = identity size `subtr` m
          zero = replicate size 0
          size = nrows m
          outUnitCircle = foldr (\(c, r) acc -> (abs c + r >= 1) || acc) False (gershgorinCircles m)  
 
doIterations m b x eps outUnitCircle cnt
    | outUnitCircle 
        && cnt == 20                 = Left 0
    | outUnitCircle
        && norm2 x' >= norm2 x + 1   = doIterations m b x' eps outUnitCircle (cnt + 1) 
    | norm2 (x `subtr` x') <= eps    = Right x'
    | otherwise                      = doIterations m b x' eps outUnitCircle 0
    where x'                         = (m `multStd` x) `add` b
 
{-leftTriangle :: (Floating a, Ord a) => Matrix a -> Matrix a
leftTriangle m = 
 
gaussZeidel :: (Floating a, Ord a) => [[a] -> [[a]] -> a -> Either Int (Matrix a)
gaussZeidel m b eps = gaussZeidel' (fromLists m) (fromLists b) eps
 
gaussZeidel' :: (Floating a, Ord a) => Matrix a -> Matrix a -> a -> Either Int (Matrix a)
gaussZeidel' m b eps = doGZ -}