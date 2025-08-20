using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()

@info("Installation done. You can now run the Atitlan.jl scripts.")
