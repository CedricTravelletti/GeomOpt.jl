using AtomsBase
using AtomsCalculators
using EmpiricalPotentials
using ExtXYZ
using Unitful
using UnitfulAtomic
using GeomOpt
using Optim

fname = joinpath(pkgdir(EmpiricalPotentials), "/home/cedric/PostPhD/Dev/DFTK/EmpiricalPotentials.jl/data", "TiAl-1024.xyz")
data = ExtXYZ.load(fname) |> FastSystem

lj = LennardJones(-1.0u"meV", 3.1u"Å",  13, 13, 6.0u"Å")

# Convert to AbstractSystem, so we have the `particles` attribute.
particles = map(data) do atom
    Atom(; pairs(atom)...)
    end
system = AbstractSystem(data; particles)

method = Optim.LBFGS()

optim_options = Optim.Options(f_tol=1e-8, g_tol=1e-8, iterations=100, show_trace=true,extended_trace=true)

results = optimize_geometry(system, lj;
            method=method, optim_options=optim_options)
println(results)
