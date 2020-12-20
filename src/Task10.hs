module Task10 where

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
import Task9
import Task3
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

multGivens :: Num a => [(Int, Int, a, a)] -> [[a]] -> [[a]]
multGivens gs m = foldr (\(i, j, c, s) acc -> givensRotation acc i j c s) m gs

{- 
QR-algorithm of calculation the matrix spectrum for the tridiagonal matrices. 
Complexity O(n^2) per iteration. Input: A, epsilon. 
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