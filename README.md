# Atitlan.jl
Zircon Data analysis of Atitlan (and more) using [MagmaThermoKinematics.jl](https://github.com/boriskaus/MagmaThermoKinematics.jl)

## Setup
To set up the julia environment and installing all the required pacakges one can run the following commands in julia:
```julia
using Pkg
julia> ]
(@v1.11) pkg> activate .
(@Atitlan) pkg> instantiate
```
or alternatively, you can use the install.jl script:
```julia
julia> include("scripts/Installation.jl")
```

## Scripts
- `AtitlanSetup.jl`: Sets up the Atitlan region, including topography and initial conditions.
- `AtitlanSetup_MTK.jl`: Sets up the Atitlan region for use with MagmaThermalKinetics.jl, plotting routines, recharge rates via multiple dispatch or the MTK routines.
- `Interpret_ZirconData_Atitlan.jl`: Interprets the model results and computes zircon ages, cumulative PDF, and other statistics useful for the analysis. Can be used independently but it is designed to run at the end of each simulation done via `AtitlanSetup_MTK.jl`.
