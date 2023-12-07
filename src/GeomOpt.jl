module GeomOpt
using LinearAlgebra
using StaticArrays
using Optim, LineSearches
using AtomsBase
using AtomsCalculators
using Unitful
using UnitfulAtomic

include("atomsbase_interface.jl")
include("optimization.jl")

end
