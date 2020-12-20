module Task12 where

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
import Task11
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


{- Test graphs on nonisomorphism. Complexity O(Time(qrEVShifts)).
    Input: A(G_1), A(G_2). -}
isomorphic :: (Floating a, Ord a) => [[a]] -> [[a]] -> Bool
isomorphic g1 g2
    | length g1 /= length g2 = False
    | otherwise              = norm2 (g1Spec `subtr` g2Spec) <= 10 * eps
    where
        eps    = 0.00001
        g1'    = toLists $ fst $ getTridiagonal g1
        g2'    = toLists $ fst $ getTridiagonal g2
        g1Spec = fromLists $ map (: []) $ DL.sort $ fst $ qrEVShifts g1' eps
        g2Spec = fromLists $ map (: []) $ DL.sort $ fst $ qrEVShifts g2' eps