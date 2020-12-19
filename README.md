# computational-linear-algebra
### Some computational linear algebra algorithms.

#### Haskell and Cabal installation commands:
``` curl --proto '=https' --tlsv1.2 -sSf https://get-ghcup.haskell.org | sh
    ghcup install ghc
    ghcup install cabal
    cabal install cabal-install
```

#### To use the module come in the root folder at first. Then enter: 
```cabal configure
   cabal build
   cabal install
```
#### Now you are ready to use the module. Enter:
```cabal repl
```
   
#### Some examples of input:
```qrEVShifts [[1, 2, 3], [1, 2, 3], [1, 2, 3]] 0.00001
   getTridiagonal [[1, 2, 3], [1, 2, 3], [1, 2, 3]]
   simpleIterationMaxEV [[1, 2, 3], [1, 2, 3], [1, 2, 3]] [[1], [1], [1]] 0.00001 100
```

#### More examples you can find in the Linalg_test.hs file.
