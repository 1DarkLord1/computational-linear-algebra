module Task8 where

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
import Task4
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

{- QR-algorithm of calculation the matrix spectrum. Complexity O(n^3) per iteration.
    Input: A, epsilon. -}
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
