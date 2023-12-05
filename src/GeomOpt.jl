module GeomOpt
using LinearAlgebra
using StaticArrays
using Optim
using AtomsBase
import AtomsCalculators
using Unitful
using UnitfulAtomic

include("atomsbase_interface.jl")
include("optimization.jl")

end
