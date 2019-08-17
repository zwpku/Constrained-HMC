# Constrained-HMC
Monte Carlo sampler on submanifolds

## Prepare
#### 1. Install the following dependent Julia libraries. 

- [HomotopyContinuation](https://www.juliahomotopycontinuation.org/), 
 [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl), used to solve multiple solutions of polynomial contraint equations when the submanifold is algebraic.
- [PolynomialRoots](https://github.com/giordano/PolynomialRoots.jl), used to
  find all root of a (scalar) polynomial equation.
- [IterativeSolvers](https://github.com/JuliaMath/IterativeSolvers.jl), used
  to solve matrix systems for veolocity constraints.
- [YAML](https://github.com/BioJulia/YAML.jl), used to parse the configure file [cfg.yml](./cfg.yml).
- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/index.html)
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

## Usage

#### 1. Enter the directory containing source files and create a working directory to run the code.

```
  	cd ./Constrained-HMC
	cp -r working_dir_template working_dir_task1
```
In the example above, the directory ./working_dir_task1 is created by duplicating the template directory [./working_dir_template](./working_dir_template)

#### 2. Create a model file in the working directory. 

It is necessary to specify the model (e.g., parameters, potential, and the mapping whose zero level-set defines the submanifold) in order to run the code. Several model files are provided under [./model_files](./model_files) directory and they can be used as templates. For example, to sample the 3D torus, we use
```
    cd working_dir_task1
    cp ../model_files/3d_torus.jl .
```

#### 3. Create the configure file cfg.yml in the working directory. 
The file cfg.yml should be prepared under the working directory. A template configure file 
[cfg.yml](./cfg.yml) is provided under [./working_dir_template](./working_dir_template).

#### 4. Set and check the configure parameters in the configure file.

In particular, the parameter "model_file_name" should be set to 
the name of the model file created in step 2. For example, the following line in ./cfg.yml tells the code that we want to use the model provided in the file ./3d_torus.jl.
      
```
      model_file_name : 3d_torus.jl
```

#### 5. Run the source code.

```
    julia ../src/constrained_hmc.jl
```

