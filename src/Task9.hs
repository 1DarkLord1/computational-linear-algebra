module Task9 where

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
import Task5
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

multHouseholderRight :: Num a => [[a]] -> [[a]] -> [[a]]
multHouseholderRight m' v' = toLists $ m `subtr` prod
    where
        m    = fromLists m'
        v    = fromLists v'
        vt   = transpose v
        prod = scaleMatrix 2 ((m `multStd` v) `multStd` vt)

{- Matrix transformation to tridiagonal form. Complexity O(n^3).
    Input: A. -}
getTridiagonal :: (Floating a, Ord a) => [[a]] -> (Matrix a, Matrix a)
getTridiagonal m = appToPair fromLists $ second trans tridiag
    where
        idm                = idMatrix $ length m
        tridiag            = foldl handler (m, idm) [1..(length m - 1)]
        handler (a, q) k 
            | norm2 v' == 0 || u == e1 = (a, q)
            | otherwise                = (newa, newq)
            where
                col  = trans a !! (k - 1)
                v'   = fromLists $ zipWith (\j x -> if j <= k then [0] else [x]) [1..length col] col
                u    = scaleMatrix (1 / norm2 v') v'
                e1   = mapPos (\(j, _) _-> if j == k + 1 then 1 else 0) v'
                v    = toLists $ scaleMatrix (1 / norm2 (u `subtr` e1)) (u `subtr` e1)
                newa = multHouseholderRight (multHouseholder a v) v
                newq = multHouseholder q v

cursorp :: ZP.Zipper a -> Int
cursorp (ZP.Zip l _) = length l

rightn :: Int -> ZP.Zipper a -> ZP.Zipper a
rightn 0 z = z
rightn n z = rightn (n - 1) (ZP.right z) 

givensRotationZ :: Num a => [ZP.Zipper a] -> Int -> Int -> a -> a -> [ZP.Zipper a]
givensRotationZ m i j c s = zipWith subst [1..] m
    where
        ui      = m !! (i - 1)
        uj      = m !! (j - 1)
        curspUi = cursorp ui
        curspUj = cursorp uj   
        uinew   = rightn curspUi $ ZP.fromList $ 
                   zipWith (\xi xj -> c * xi + s * xj) (ZP.toList ui) (ZP.toList uj)
        ujnew   = rightn curspUj $ ZP.fromList $ 
                   zipWith (\xi xj -> (- s) * xi + c * xj) (ZP.toList ui) (ZP.toList uj)
        subst pos row
            | pos == i  = uinew
            | pos == j  = ujnew
            | otherwise = row

{- QR decomposition for the tridiagonal matrices. Complexity O(n^2).
    Input: A, epsilon. -}
qrDecompTridiagonal :: (Floating a, Ord a) => [[a]] -> ([(Int, Int, a, a)], Matrix a)
qrDecompTridiagonal m = second (fromLists . fmap ZP.toList) $ 
                        foldl handler ([], zipped_m) [1..(length m - 1)]
    where
        zipped_m = fmap ZP.fromList m
        handler (gs, r) k
            | n == 0    = (gs, r)
            | otherwise = ((k + 1, k, c, s) : gs, map ZP.right $ givensRotationZ r (k + 1) k c s)
            where
                col = map ZP.cursor r
                xi  = col !! (k - 1)
                xj  = col !! k
                n   = sqrt $ xi * xi + xj * xj
                c   = xi / n
                s   = (- xj) / n