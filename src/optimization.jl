# 
# Convenience functions for computation of forces and energies.
#
export energy_wrt_positions, forces_wrt_positions, construct_optimization_problem

abstract type AbstractCalculator end

# TODO: deprecated.
function energy_wrt_positions(system::AbstractSystem, positions::AbstractVector{<:AbstractVector{T}}, calculator::AbstractCalculator) where {T <: Real}
    new_system = update_positions(system, positions)
    AtomsCalculators.potential_energy(new_system, calculator)
end

# TODO: deprecated.
function forces_wrt_positions(system::AbstractSystem, positions::AbstractVector{<:AbstractVector{T}}, calculator::AbstractCalculator) where {T <: Real}
    new_system = update_positions(system, positions)
    AtomsCalculators.forces(new_system, calculator)
end

function construct_optimization_problem(system::AbstractSystem, calculator::AbstractCalculator; cartesian=false, units=u"angstrom")
    f = if cartesian
        function(x::AbstractVector{<:Real})
            # Attach unit information.
            x = x .* units
            new_system = update_optimizable_coordinates_cart(system, x)
            AtomsCalculators.potential_energy(new_system, calculator)
        end
    else
        function(x::AbstractVector{<:Real})
            new_system = update_optimizable_coordinates(system, x)
            AtomsCalculators.potential_energy(new_system, calculator)
        end
    end
    return f
end

function fg!(F,G,x)
    energy = AtomsCalculators.potential_energy(new_system, calculator)
  if G != nothing
      forces = AtomsCalculators.forces(system, calculator; kwargs..., precomputed=true)
      G .= forces
  end
  if F != nothing
    return energy
  end
end
