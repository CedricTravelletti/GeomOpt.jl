using DFTK: AtomsCalculators
#= Test Geometry Optimization on an aluminium supercell.
=#
using Printf
using LinearAlgebra
using DFTK
using ASEconvert
using Unitful
using UnitfulAtomic
using Random
using Optim, LineSearches

using GeomOpt


function build_al_supercell(rep=1)
    a = 7.65339 # true lattice constant.
    lattice = a * Matrix(I, 3, 3)
    Al = ElementPsp(:Al; psp=load_psp("hgh/lda/al-q3"))
    atoms     = [Al, Al, Al, Al]
    positions = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
    unit_cell = periodic_system(lattice, atoms, positions)

    # Make supercell in ASE:
    # We convert our lattice to the conventions used in ASE, make the supercell
    # and then convert back ...
    supercell_ase = convert_ase(unit_cell) * pytuple((rep, 1, 1))
    supercell     = pyconvert(AbstractSystem, supercell_ase)

    # Unfortunately right now the conversion to ASE drops the pseudopotential information,
    # so we need to reattach it:
    supercell = attach_psp(supercell; Al="hgh/lda/al-q3")
    return supercell
end;

al_supercell = build_al_supercell(1)

# Set everything to optimizable.
al_supercell = clamp_atoms(al_supercell, nothing)

# Create a simple calculator for the model.
calculator = DFTKCalculator(;Ecut=30.0, kgrid=[4, 4, 4], tol=1e-6, temperature=1e-4, verbose=true)

# Starting position is a random perturbation of the equilibrium one.
Random.seed!(1234)
x0 = vcat(get_positions_cart(al_supercell)...)
σ = 0.3u"angstrom"; x0_pert = x0 + σ * rand(Float64, size(x0))

energy_true = AtomsCalculators.potential_energy(al_supercell, calculator)
energy_pert = AtomsCalculators.potential_energy(
                update_optimizable_coordinates_cart(al_supercell, x0_pert), calculator)
@printf "Initial guess distance (norm) from true parameters %.3e bohrs.\n" austrip(norm(x0 - x0_pert))
@printf "Initial regret %.3e.\n" energy_pert - energy_true

# linesearch = LineSearches.BackTracking(order=2, ρ_lo=0.01, ρ_hi=0.5, c_1=100.0, maxstep=1.0, iterations=2)
linesearch = LineSearches.Static()
method = Optim.BFGS(;alphaguess=InitialStatic(alpha=10), linesearch=linesearch)
results = optimize_geometry(al_supercell, calculator, austrip.(x0_pert);
            method=method,
            optim_options=Optim.Options(f_tol=1e-32, iterations=10, show_trace=true,extended_trace=true))
println(results)
