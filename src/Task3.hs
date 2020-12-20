module Task3 where

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

{- Givens rotation of matrix. Complexity O(n).
    Input: A, i, j, c, s. -}
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
