# Constrained-HMC
Monte Carlo sampler on submanifolds

## PREPARE
### 1. Install dependent Julia libraries, if necessary. 

The following Julia packages are required.

- [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/)
- DynamicPolynomials 
- PolynomialRoots 
- LinearAlgebra 
- IterativeSolvers 
- PyPlot

For example, HomotopyContinuation.jl package can be installed as follows.

```
  julia> import Pkg
  julia> Pkg.add("HomotopyContinuation")
```

2. Download the source code to a local directory.

```
	git clone https://github.com/zwpku/Constrained-HMC.git
```

   The code should be available in the directory ./Constrained-HMC


## USAGE

1. Enter the directory containing source files. 

```
  	cd ./Constrained-HMC
```

