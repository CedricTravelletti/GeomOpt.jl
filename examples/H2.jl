#= Test Geometry Optimization on an aluminium supercell.
=#
using Printf
using LinearAlgebra
using DFTK
using Unitful
using UnitfulAtomic
using Optim

using GeomOpt


kgrid = [1, 1, 1]       # k-point grid
Ecut = 5.0                # kinetic energy cutoff in Hartree
tol = 1e-8              # tolerance for the optimization routine
a = 10.                  # Big box around the atoms.
lattice = a * I(3)
H = ElementPsp(:H; psp=load_psp("hgh/lda/h-q1"));
atoms = [H, H];
positions = [[0, 0, 0], [0, 0, .15]]
system = periodic_system(lattice, atoms, positions)

# Set everything to optimizable.
system = clamp_atoms(system, nothing)

# Create a simple calculator for the model.
calculator = DFTKCalculator(system; Ecut=Ecut, kgrid=kgrid, tol=tol, verbose=true)

method = Optim.NelderMead()
optim_options = Optim.Options(f_tol=1e-32, iterations=20, show_trace=true,extended_trace=true)

results = optimize_geometry(system, calculator; method=method, optim_options=optim_options)
println(results)
@printf "Bond length: %3f bohrs.\n" norm(results.minimizer[1:3] - results.minimizer[4:end])
