module Task1 where

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


{- Simple iteration method of solving linear equations systems. Complexity O(n^2) per iteration.
    Input: A, b, epsilon.
 -}
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