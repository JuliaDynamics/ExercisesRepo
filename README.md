**Contents**
- [Tutorials](#tutorials)
- [Installation and reproducibility](#installation-and-reproducibility)
- [Exercise datasets](#exercise-datasets)
- [Multiple choice questions](#multiple-choice-questions)
- [Figures](#figures)
- [Interactive applications, videos, extra plots](#interactive-applications--videos--extra-plots)

## Introduction
Here some introductory text about this repository will be stated.

## Tutorials for Julia and related packages

- https://www.youtube.com/watch?v=8h8rQyEpiZA&t=42s : A beginner introduction to the Julia language
- https://www.youtube.com/watch?v=Fi7Pf2NveH0 : Intensive Julia workshop, for people already familiar with programming that want to transition to Julia
- https://www.youtube.com/watch?v=13hqE_1a158 : Introduction to DynamicalSystems.jl
- https://juliadynamics.github.io/DynamicalSystems.jl/dev/ : Documentation for DynamicalSystems.jl
- https://juliadynamics.github.io/JuliaDynamics/ : Website of the JuliaDynamics organization
- https://www.youtube.com/watch?v=KPEqYtEd-zY : Introduction to solving differential equations in Julia, which is generally useful for nonlinear dynamics
- https://julialang.org/community/ : Resources for asking questions about Julia
- https://discourse.julialang.org/ : Official Julia forum (and also the main platform that newcomers ask questions)


## Installation and reproducibility

The accompanying code base used here is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ExercisesRepo

To (locally) reproduce this project, first install Julia and then do the following:

0. Download this repository as-is and export it to some folder.
0. Install `DrWatson` on your general Julia installation by doing:
   ```
   julia> using Pkg; Pkg.add("DrWatson")
   ```
1. Then do:
   ```
   julia> Pkg.activate("path/to/the/downloaded/project")
   julia> Pkg.instantiate()
   ```
1. This repository uses the Python library matplotlib (package `PyPlot`) for plots. You could use any plotting package instead, but if you want to use `PyPlot` and replicate exactly our plots, then run the following commands to ensure a working installation for all operating systems:
   ```
   julia> ENV["PYTHON"] = ""
   julia> Pkg.add("PyCall"); Pkg.build("PyCall")
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
The datasets that are used in the book exercises are contained in the `exercise_data` folder, all being in the same text-based format. The same folder contains information of where this data is coming from: `data_explanations.md`. Some data are generated from simulations in the script `generating_code.jl`.

## Multiple choice questions
Multiple choice questions that we use during lecturing to increase student involvement are in the `multiple_choice` folder.

## Figures
The code that creates the figures of our book is in the `figure_generation` folder. Notice however that some figures were made with (or enhanced by) PowerPoint and thus we do not share this here.

## Interactive applications, videos, extra plots
In the folder `animations_code` we provide animation-related scripts (some of which are self-documented by comments). These provide the following extras that we use while lecturing the various chapters:
* launch interactive applications
* create videos

We have also recorded videos for easier access, stored in `animations_videos`.
