module Task6 where

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
import Task5
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

{- Householder reflection method of QR decomposition. Complexity O(n^3). -}
qrDecompHouseholder :: (Floating a, Ord a) => [[a]] -> (Matrix a, Matrix a)
qrDecompHouseholder m = appToPair fromLists $ first trans qr
    where
        idm                = idMatrix $ length m
        qr                 = foldl handler (idm, m) [1..(length m - 1)]
        handler (q, r) k 
            | norm2 v' == 0 || u == e1 = (q, r)
            | otherwise                = (multHouseholder q v, multHouseholder r v)
            where
                col = trans r !! (k - 1)
                v'  = fromLists $ zipWith (\j x -> if j < k then [0] else [x]) [1..length col] col
                u   = scaleMatrix (1 / norm2 v') v'
                e1  = mapPos (\(j, _) _-> if j == k then 1 else 0) v'
                v   = toLists $ scaleMatrix (1 / norm2 (u `subtr` e1)) (u `subtr` e1)
 