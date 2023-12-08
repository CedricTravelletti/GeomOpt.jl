#=
For the moment, optimization is delegated to Optim.jl (hardcoded dependency). 
In order for Optim.jl to use gradients that are returned directly upon function 
evaluation, one has to use the fix described at: 
https://github.com/JuliaNLSolvers/Optim.jl/blob/master/docs/src/user/tipsandtricks.md

This is the reason for the cumbersome implementation of optimization with gradients.
=#

export construct_optimization_function, construct_optimization_function_w_gradients, optimize_geometry


""" 
    By default we work in cartesian coordinaes.
    Note that internally, when optimizing the cartesian positions, atomic units 
    are used.

"""
function construct_optimization_function(system::AbstractSystem, calculator::AbstractCalculator; cartesian=true, kwargs...)
    f = if cartesian
        function(x::AbstractVector{<:Real})
            x = 1u"bohr" .* x # Work in atomic units.
            new_system = update_optimizable_coordinates_cart(system, x)
            AtomsCalculators.potential_energy(new_system, calculator, calculator.state; kwargs...)
        end
    else
        function(x::AbstractVector{<:Real})
            new_system = update_optimizable_coordinates(system, x)
            AtomsCalculators.potential_energy(new_system, calculator, calculator.state; kwargs...)
        end
    end
    return f
end

function construct_optimization_function_w_gradients(system::AbstractSystem, calculator::AbstractCalculator; cartesian=true, kwargs...)
    fg! = function(F::Union{Nothing, Real}, G::Union{Nothing, AbstractVector{<:Real}}, x::AbstractVector{<:Real})
        x = 1u"bohr" .* x # Work in atomic units.
        new_system = update_optimizable_coordinates_cart(system, x)
        energy = AtomsCalculators.potential_energy(new_system, calculator, calculator.state; kwargs...)
        if G != nothing
            forces = AtomsCalculators.forces(new_system, calculator, calculator.state; kwargs...)
            # Translate the forces vectors on each particle to a single gradient for the optimization parameter.
            forces_concat = mask_vector_list(forces, get_optimizable_mask(new_system))
            G .= forces_concat
        end
        if F != nothing
          return energy
        end
    end
    return fg!
end

function optimize_geometry(system::AbstractSystem, calculator::AbstractCalculator, x0::AbstractVector{<:Real};
        no_gradients=false,
        method=Optim.NelderMead(),
        optim_options=Optim.Options(show_trace=true,extended_trace=true), kwargs...)
    x0 = Vector(x0) # Optim modifies x0 in-place, so need a mutable type.

    if no_gradients
        f = construct_optimization_function(system, calculator; kwargs...)
    else
        fg! = construct_optimization_function_w_gradients(system, calculator; kwargs...)
        f = Optim.only_fg!(fg!)
    end
    optimize(f, x0, method, optim_options; kwargs...)
end
function optimize_geometry(system::AbstractSystem, calculator::AbstractCalculator;
        no_gradients=false,
        method=Optim.NelderMead(),
        optim_options=Optim.Options(show_trace=true,extended_trace=true), kwargs...)
    # Use current system parameters as starting positions.
    x0 = austrip.(get_optimizable_coordinates_cart(system)) # Optim modifies x0 in-place, so need a mutable type.

    if no_gradients
        f = construct_optimization_function(system, calculator; kwargs...)
    else
        fg! = construct_optimization_function_w_gradients(system, calculator; kwargs...)
        f = Optim.only_fg!(fg!)
    end
    optimize(f, x0, method, optim_options; kwargs...)
end
