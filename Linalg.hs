{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE InstanceSigs #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE UndecidableInstances #-}

import Data.Matrix

type MatrixDouble = Matrix Double

subtr :: Num a => Matrix a -> Matrix a -> Matrix a
subtr m1 m2 = elementwise (\x y -> x - y) m1 m2

add :: Num a => Matrix a -> Matrix a -> Matrix a
add m1 m2 = elementwise (\x y -> x + y) m1 m2 

norm2 :: (Floating a, Num a) => Matrix a -> a
norm2 v = sqrt $ foldr (\x r -> x * x + r) 0 (toList v) 

simpleIteration :: [[Double]] -> [[Double]] -> Double -> Either Int MatrixDouble
simpleIteration m b eps = simpleIteration' (fromLists m) (fromLists b) eps

simpleIteration' :: MatrixDouble -> MatrixDouble -> Double -> Either Int MatrixDouble
simpleIteration' m b eps = doIterations m' b (fromList size 1 zero) eps 0 
	where 
		  m'   = identity size `subtr` m
		  zero = replicate size 0
		  size = nrows m


doIterations m b x eps cnt
	| cnt == 20 				  = Left 0
	| norm2 x' >= (norm2 x) + 1   = doIterations m b x' eps (cnt + 1) 
	| norm2 (x `subtr` x') <= eps = Right x'
	| otherwise                   = doIterations m b x' eps cnt
	where x'                      = (m `multStd` x) `add` b 
