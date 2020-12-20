module Task11 where

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
import Task10
import Task9
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
Complexity O(n^2) per iteration. Input: A, epsilon.
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
