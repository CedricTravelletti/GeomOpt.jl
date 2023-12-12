#
# Interface between AtomsBase and GeomOpt that provides 
# utility functions for manipulating systems.
#
# This interface helps defining which particles in the system are considered 
# optimizable. Note that by default all particles are assumed to be optimizable. 
# The user can call clamp_atoms to fix atoms whose position should not be optimized.
#
# IMPORTANT: Note that we always work in cartesian coordinates.
#
export fractional_to_cartesian, cartesian_to_fractional, update_positions
export update_optimizable_coordinates, set_optimizable_mask
export get_optimizable_mask, get_optimizable_coordinates, mask_vector_list, clamp_atoms

@doc raw"""
    fractional_to_cartesian(system::AbstractSystem, positions::AbstractVector{<:AbstractVector{<:Real}})

Given a list of fractional coordinates, convert them to cartesian.
"""
function fractional_to_cartesian(system::AbstractSystem, positions::AbstractVector{<:AbstractVector{<:Real}})
    [fractional_to_cartesian(system, position) for position in positions]
end
function fractional_to_cartesian(system::AbstractSystem, positions::AbstractVector{<:Real})
    # Get lattice matrix for converting fractional to cartesian.
    lattice = hcat(system.bounding_box...)
    lattice * positions
end

@doc raw"""
    cartesian_to_fractional(system::AbstractSystem, positions::AbstractVector{<:AbstractVector{<:Real}})

Given a list of cartesian coordinates, convert them to fractional.
"""
function cartesian_to_fractional(system::AbstractSystem, positions::AbstractVector{<:AbstractVector{<:Unitful.Length}})
    [cartesian_to_fractional(system, position) for position in positions]
end
function cartesian_to_fractional(system::AbstractSystem, positions::AbstractVector{<:Unitful.Length})
    # Get lattice matrix for converting fractional to cartesian.
    lattice = hcat(system.bounding_box...)
    # Make sure we are working in the same units. Have to strip for the dimensionless product.
    inv(ustrip(uconvert.(u"angstrom", lattice))) * ustrip(uconvert.(u"angstrom",positions))
end

@doc raw"""
    update_postions(system::AbstractSystem, positions::Vector{<:AbstractVector{<:Real}}) where {L <: Unitful.Length}

Creates a new system based on ``system`` but with atoms positions updated to the ones specified in `positions`, 
using cartesian coordinates.
"""
function update_positions(system::AbstractSystem, positions::AbstractVector{<:AbstractVector{L}}) where {L <: Unitful.Length}
    particles = [
                 Atom(atom, position=SVector{3,L}(position)) 
                 for (atom, position) in zip(system.particles, positions)]
    AbstractSystem(system, particles=particles)
end

@doc raw"""
    update_optimizable_coordinates(system::AbstractSystem, positions::Vector{<:AbstractVector{<:L}}) where {L <: Unitful.Length}

Creates a new system based on ``system`` with the optimizable coordinates are
updated to the ones provided (in the order in which they appear in system.particles), cartesian coordinates version..
"""
function update_optimizable_coordinates(system::AbstractSystem, positions::AbstractVector{<:L}) where {L <: Unitful.Length}
    mask = get_optimizable_mask(system)
    new_positions = Vector.(position(system)) # make mutable.
    new_positions[mask] = [eachcol(reshape(positions, 3, :))...]
    update_positions(system, new_positions)
end

@doc raw"""
    set_optimizable_mask(system::AbstractSystem, mask::AbstractVector{<:Bool})

Sets the mask defining which coordinates of the system can be optimized. 
The mask is a vector of booleans, specifying which of the atoms can be optimized.

By default (when no mask is specified), all particles are assumed optimizable.

"""
function set_optimizable_mask(system::AbstractSystem, mask::AbstractVector{<:Bool})
    particles = [Atom(atom; optimizable=m) for (atom, m) in zip(system.particles, mask)]
    AbstractSystem(system, particles=particles)
end

@doc raw"""
    get_optimizable_mask(system::AbstractSystem) -> AbstractVector{<:Bool}

Returns the optimizable mask of the system (see documentation for `set_optimizable_mask`.
"""
function get_optimizable_mask(system::AbstractSystem)
    # If flag not set, the atom is considered to be optimizable.
    [haskey(a, :optimizable) ? a[:optimizable] : true for a in system.particles]
end

@doc raw"""
    Given a list of vectors and a mask, return the masked version in a one-dimensional list.
    This is used to concatenate the forces returned in the scfres into a single list, 
    that can be used by the optimizer.

    """
function mask_vector_list(x::AbstractVector{<:AbstractVector{<:Any}}, mask::AbstractVector{<:Bool})
    collect(Iterators.flatten(x[mask]))
end

function get_optimizable_coordinates(system::AbstractSystem)
    mask = get_optimizable_mask(system)
    return collect(Iterators.flatten(system[mask, :position]))
end

@doc raw"""
    Clamp given atoms if the system. Clamped atoms are fixed and their positions 
    will not be optimized. The atoms to be clamped should be given as a list of 
    indies corresponding to their positions in system.particles.

    """
function clamp_atoms(system::AbstractSystem, clamped_indexes::Union{AbstractVector{<:Integer},Nothing})
    mask = trues(length(system.particles))
    mask[clamped_indexes] .= false
    clamped_system = set_optimizable_mask(system, mask)
    return clamped_system
end
