module Task13 where

{-# LANGUAGE RankNTypes,
             FlexibleContexts,
             FlexibleInstances,
             MultiParamTypeClasses,
             UndecidableInstances,
             AllowAmbiguousTypes,
             DataKinds,
             ScopedTypeVariables,
#-}

import Task9
import Task11
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

getAdjMap :: Num a => Int -> DM.Map Int (DM.Map Int a)
getAdjMap n = adjMap 
    where
        buildRow = DM.fromList $ map (\k -> (k, 0)) [0..n - 1]
        adjMap   = DM.fromList $ map (\k -> (k, buildRow)) [0..n - 1]

addEdge (x, y) = DM.update (Just . DM.update (Just . (+ 1)) y) x

buildAdjMatrix :: Num a => Int -> (Int -> DM.Map Int (DM.Map Int a) -> DM.Map Int (DM.Map Int a)) -> Matrix a
buildAdjMatrix n genEdges = m
    where
        adjMap      = getAdjMap n
        adjMap'     = foldr genEdges adjMap [0..n - 1]
        m           = matrix n n fill
        fill (i, j) = x
            where
                i' = i - 1
                j' = j - 1
                x  = fromJust $ DM.lookup (j - 1) $ fromJust $ DM.lookup (i - 1) adjMap'

buildGraph1 :: Int -> Matrix Double
buildGraph1 n = buildAdjMatrix (n * n) genEdges
    where
        genEdges x = foldr (\y acc -> genEdges' (x, y) . acc) id [0..n - 1]
        genEdges' (x, y)
            | x >= n    = id
            | otherwise = addEdge (v, u1) . addEdge (v, u2) . 
                          addEdge (v, u3) . addEdge (v, u4) .
                          addEdge (v, u5) . addEdge (v, u6) . 
                          addEdge (v, u7) . addEdge (v, u8)
            where
                x1 = (x + 2 * y) `mod` n
                x2 = (x - 2 * y + 3 * n) `mod` n
                x3 = (x + 2 * y + 1) `mod` n
                x4 = (x - 2 * y - 1 + 3 * n) `mod` n
                y1 = (y + 2 * x) `mod` n
                y2 = (y - 2 * x + 3 * n) `mod` n
                y3 = (y + 2 * x + 1) `mod` n
                y4 = (y - 2 * x - 1 + 3 * n) `mod` n
                u1 = x1 * n + y
                u2 = x2 * n + y
                u3 = x3 * n + y
                u4 = x4 * n + y
                u5 = x * n + y1
                u6 = x * n + y2
                u7 = x * n + y3
                u8 = x * n + y4
                v  = x * n + y
                
inv :: Int -> Int -> Int
inv x md = ((xinv `mod` md) + md) `mod` md
    where 
        xinv = snd $ gcdExt x md  

buildGraph2 :: Int -> Matrix Double
buildGraph2 p = buildAdjMatrix (p + 1) genEdges
    where
        genEdges x  
            | x == p    = addEdge (x, xInv) . addEdge (x, x) . addEdge (x, x)
            | otherwise = addEdge (x, xInv) . addEdge (x, x1) . addEdge (x, x2)
            where
                xInv
                    | x == 0    = p
                    | x == p    = 0
                    | otherwise = inv x p
                x1 = (x + 1) `mod` p
                x2 = (x - 1 + p) `mod` p

{- Calculation the optimal alpha for expander. Complexity O(Time(qrEVShifts)). -}
expanderAlpha :: Int -> Matrix Double -> Double -> Double
expanderAlpha n g d = max (abs ev1) (abs ev2) / d 
    where
        eps       = 0.00001
        evs       = fst $ qrEVShifts (toLists $ fst $ getTridiagonal $ toLists g) eps
        sortedEVs = reverse $ DL.sort evs
        ev1       = head $ tail sortedEVs
        ev2       = last sortedEVs

{- Optimal alpha for the first graph. Input: n. -}
expanderAlpha1 :: Int -> Double 
expanderAlpha1 n = expanderAlpha n g 8
    where
        g = buildGraph1 n

{- Optimal alpha for the second graph. Input: p. -}
expanderAlpha2 :: Int -> Double
expanderAlpha2 n = expanderAlpha (n + 1) g 3
    where
        g = buildGraph2 n