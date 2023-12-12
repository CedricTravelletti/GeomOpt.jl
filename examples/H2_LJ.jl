using EmpiricalPotentials: AtomsCalculators
#= Test Geometry Optimization on an aluminium supercell.
=#
using Printf
using LinearAlgebra
using EmpiricalPotentials
using Unitful
using UnitfulAtomic
using Optim

using GeomOpt


a = 10.                  # Big box around the atoms.
lattice = a * I(3)
H = ElementPsp(:H; psp=load_psp("hgh/lda/h-q1"));
atoms = [H, H];
positions = [[0, 0, 0], [0, 0, .19]]
system = periodic_system(lattice, atoms, positions)

lj = LennardJones(-1.17u"hartree", 0.743u"angstrom", 1, 1, 0.6u"nm")

method = Optim.LBFGS()
optim_options = Optim.Options(f_tol=1e-16, iterations=200, show_trace=true,extended_trace=true)

results = optimize_geometry(system, lj; method=method, optim_options=optim_options)
println(results)
@printf "Bond length: %3f bohrs.\n" norm(results.minimizer[1:3] - results.minimizer[4:end])
