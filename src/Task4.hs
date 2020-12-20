module Task4 where

{-# LANGUAGE RankNTypes,
             FlexibleContexts,
             FlexibleInstances,
             MultiParamTypeClasses,
             UndecidableInstances,
             AllowAmbiguousTypes,
             DataKinds,
             ScopedTypeVariables,
#-}


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

trans :: [[a]] -> [[a]]
trans = toLists . transpose . fromLists

idMatrix :: Num a => Int -> [[a]]
idMatrix sz = toLists $ identity sz

firstNonzero :: (Num a, Ord a) => [a] -> Int -> Int
firstNonzero v k = foldl (\acc (j, x) -> if j >= k && x /= 0 && acc == 0 
                                       then j 
                                       else acc) 0 (zip [1..] v) 

appToPair f (x, y) = (f x, f y)

{- Givens rotation method of QR decomposition. Complexity O(n^3).
    Input: A. -}
qrDecompGivens :: (Floating a, Ord a) => [[a]] -> (Matrix a, Matrix a)
qrDecompGivens m = appToPair fromLists $ first trans qr
    where
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
                    | j <= i    = ((q', r'), xi)
                    | otherwise = ((givensRotation q' j i c s, givensRotation r' j i c s), (- s) * xj + c * xi)
                        where
                            n = sqrt $ xi * xi + xj * xj
                            c = xi / n
                            s = (- xj) / n