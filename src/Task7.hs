module Task7 where

{-# LANGUAGE RankNTypes,
             FlexibleContexts,
             FlexibleInstances,
             MultiParamTypeClasses,
             UndecidableInstances,
             AllowAmbiguousTypes,
             DataKinds,
             ScopedTypeVariables,
#-}

import Task1
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

{- Simple iteration method of calculation the maximum modulo eigenvalue. Complexity O(n^2) per iteration.
    Input: A, x0, epsilon, maximal allowed count of iterations. -}
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