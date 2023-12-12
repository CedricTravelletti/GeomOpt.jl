#=
For the moment, optimization is delegated to Optim.jl (hardcoded dependency). 
In order for Optim.jl to use gradients that are returned directly upon function 
evaluation, one has to use the fix described at: 
https://github.com/JuliaNLSolvers/Optim.jl/blob/master/docs/src/user/tipsandtricks.md

This is the reason for the cumbersome implementation of optimization with gradients.

Note that by default all particles in the system are assumed optimizable.

IMPORTANT: Note that we always work in cartesian coordinates.

=#

# TODO: Optimization with only_fg! does not seem to be working (function is re-called upon gradient computation). 
# Given that the initial guess in the SCF now reduces that computation time to nothing, do we really want 
# to use this convoluted syntax?

export construct_optimization_function, construct_optimization_function_w_gradients, optimize_geometry


""" 
    By default we work in cartesian coordinaes.
    Note that internally, when optimizing the cartesian positions, atomic units 
    are used.

"""
function construct_optimization_function(system::AbstractSystem, calculator; kwargs...)
    f = function(x::AbstractVector{<:Real})
            x = 1u"bohr" .* x # Work in atomic units.
            new_system = update_optimizable_coordinates(system, x)
            austrip(AtomsCalculators.potential_energy(new_system, calculator; kwargs...))
    end
    return f
end

function construct_optimization_function_w_gradients(system::AbstractSystem, calculator; kwargs...)
    fg! = function(F::Union{Nothing, Real}, G::Union{Nothing, AbstractVector{<:Real}}, x::AbstractVector{<:Real})
        x = 1u"bohr" .* x # Work in atomic units.
        new_system = update_optimizable_coordinates(system, x)
        energy = AtomsCalculators.potential_energy(new_system, calculator; kwargs...)
        if G != nothing
            forces = AtomsCalculators.forces(new_system, calculator; kwargs...)
            # Translate the forces vectors on each particle to a single gradient for the optimization parameter.
            forces_concat = mask_vector_list(forces, get_optimizable_mask(new_system))
            # NOTE: minus sign since forces are opposite to gradient.
            G .= - austrip.(forces_concat)
        end
        if F != nothing
            return austrip(energy)
        end
    end
    return fg!
end

function optimize_geometry(system::AbstractSystem, calculator, x0::AbstractVector{<:Real};
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
function optimize_geometry(system::AbstractSystem, calculator;
        no_gradients=false,
        method=Optim.NelderMead(),
        optim_options=Optim.Options(show_trace=true,extended_trace=true), kwargs...)
    # Use current system parameters as starting positions.
    x0 = Vector(austrip.(get_optimizable_coordinates(system))) # Optim modifies x0 in-place, so need a mutable type.

    if no_gradients
        f = construct_optimization_function(system, calculator; kwargs...)
    else
        fg! = construct_optimization_function_w_gradients(system, calculator; kwargs...)
        f = Optim.only_fg!(fg!)
    end
    optimize(f, x0, method, optim_options; kwargs...)
end
