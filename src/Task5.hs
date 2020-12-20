module Task5 where

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

{- Householder matrix multiplication. Complexity O(n^2).
    Input: A, v. -}
multHouseholder :: Num a => [[a]] -> [[a]] -> [[a]]
multHouseholder m' v' = toLists $ m `subtr` prod
    where
        m    = fromLists m'
        v    = fromLists v'
        vt   = transpose v
        prod = scaleMatrix 2 v `multStd` (vt `multStd` m)