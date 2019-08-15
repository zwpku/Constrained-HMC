# Constrained-HMC
Monte Carlo sampler on submanifolds

## PREPARE
#### 1. Install the following dependent Julia libraries. 

- [HomotopyContinuation](https://www.juliahomotopycontinuation.org/)
- [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl)
- [PolynomialRoots](https://github.com/giordano/PolynomialRoots.jl)
- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/index.html)
- [IterativeSolvers](https://github.com/JuliaMath/IterativeSolvers.jl)
- [PyPlot](https://github.com/JuliaPy/PyPlot.jl)

For example, the HomotopyContinuation package can be installed as follows.

```
  julia> import Pkg
  julia> Pkg.add("HomotopyContinuation")
```

#### 2. Download the source code to a local directory.

```
	git clone https://github.com/zwpku/Constrained-HMC.git
```

   Afterwards, the code should be available in the directory ./Constrained-HMC

## USAGE

1. Enter the directory containing source files. 

```
  	cd ./Constrained-HMC
```


