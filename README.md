Some info lalala.

**Contents**
- [Installation and reproducibility](#installation-and-reproducibility)
- [Exercise datasets](#exercise-datasets)

## Installation and reproducibility

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ExercisesRepo

To (locally) reproduce this project, first install Julia and then do the following:

0. Download this code base as-is.
0. Install `DrWatson` on your general Julia installation by doing:
   ```
   julia> using Pkg; Pkg.add("DrWatson")
   ```
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```
1. This repository uses matplotlib (package `PyPlot`) for plots. You could use any plotting package instead, but if you want to use `PyPlot` run the following commands to ensure a working installation for all operating systems:
   ```
   julia> ENV["PYTHON"] = ""
   julia> Pkg.add("PyCall"); Pkg.build("Pycall")
   julia> Pkg.add("PyPlot"); using PyPlot
   ```

Now all necessary packages are installed and all scripts should run out of the box.
As you will notice, all scripts start with the commands:
```julia
using DrWatson
@quickactivate "ExercisesRepo"
```
which ensures that only local directories will be used, as well as the *exact* package versions contained within the repository, leading to full reproducibility.

## Exercise datasets
The datasets that are used in the book exercises are contained in the `timeseries` folder, all being in the same text-based format. The same folder contains information of where this data is coming from: `data_explanations.md`. Some data are generated from simulations in the script `generating_code.jl`.

## Multiple choice questions
