#
# Interface between AtomsBase and GeomOpt that provides 
# utility functions for manipulating systems.
#
export fractional_to_cartesian, cartesian_to_fractional, update_positions, update_positions_cart
export update_optimizable_coordinates_cart, update_optimizable_coordinates, set_optimizable_mask, get_positions, get_positions_cart
export n_optimizable_coordinates, get_optimizable_mask, get_optimizable_coordinates_cart, mask_vector_list, clamp_atoms

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
    update_postions_cart(system::AbstractSystem, positions::Vector{<:AbstractVector{<:Real}}) where {L <: Unitful.Length}

Creates a new system based on ``system`` but with atoms positions updated to the ones specified in `positions`, 
using cartesian coordinates.
"""
function update_positions_cart(system::AbstractSystem, positions::AbstractVector{<:AbstractVector{L}}) where {L <: Unitful.Length}
    particles = [
                 Atom(atom, position=SVector{3,L}(position)) 
                 for (atom, position) in zip(system.particles, positions)]
    AbstractSystem(system, particles=particles)
end

@doc raw"""
    update_postions(system::AbstractSystem, positions::Vector{<:AbstractVector{<:Real}})

Creates a new system based on ``system`` but with atoms positions updated to the ones specified in `positions`.
Note that we work in fractional coordinates.
"""
function update_positions(system::AbstractSystem, positions::AbstractVector{<:AbstractVector{<:Real}})
    update_positions_cart(system, fractional_to_cartesian(system, positions))
end


@doc raw"""
    update_optimizablecoordinates_cart(system::AbstractSystem, positions::Vector{<:AbstractVector{<:L}}) where {L <: Unitful.Length}

Creates a new system based on ``system`` with the optimizable coordinates are
updated to the ones provided (in the order in which they appear in system.particles), cartesian coordinates version..
"""
function update_optimizable_coordinates_cart(system::AbstractSystem, positions::AbstractVector{<:L}) where {L <: Unitful.Length}
    new_positions = Vector{Vector{L}}() # Empty array of atoms for new system.
    positions = Vector(positions) # Make it mutable to use pop!.
    for atom in system.particles
        # Check if the `optimizable`key exists.
        # If yes, update positions that can be updated.
        if haskey(atom, :optimizable)
            new_position = Vector{L}()
            for k in eachindex(atom[:optimizable])
                if atom[:optimizable][k]
                    # Can be optimized, so update position.
                    push!(new_position, popfirst!(positions))
                else
                    push!(new_position, atom.position[k])# Cannot be optimized, so use current position.
                end
            end
        else
            # otherwise re-use old position.
            new_position = atom.position
        end
        push!(new_positions, new_position)
    end
    # Safety check.
    if !isempty(positions)
        error("Size of new positions and number of optimizable positions do not match")
    end
    update_positions_cart(system, new_positions)
end

@doc raw"""
    update_optimizable_coordinates(system::AbstractSystem, positions::Vector{<:AbstractVector{L}}) where {L <: Real}

Creates a new system based on ``system`` with the optimizable coordinates are
updated to the ones provided (in the order in which they appear in system.particles).
"""
function update_optimizable_coordinates(system::AbstractSystem, positions::AbstractVector{L}) where {L <: Real}
    new_positions = Vector{Vector{L}}() # Empty array of atoms for new system.
    positions = Vector(positions) # Make it mutable to use pop!.
    for atom in system.particles
        # Check if the `optimizable`key exists.
        # If yes, update positions that can be updated.
        if haskey(atom, :optimizable)
            new_position = Vector{L}()
            for k in eachindex(atom[:optimizable])
                if atom[:optimizable][k]
                    # Can be optimized, so update position.
                    push!(new_position, popfirst!(positions))
                else
                    push!(new_position, atom.position[k])# Cannot be optimized, so use current position.
                end
            end
        else
            # otherwise re-use old position.
            new_position = atom.position
        end
        push!(new_positions, new_position)
    end
    # Safety check.
    if !isempty(positions)
        error("Size of new positions and number of optimizable positions do not match")
    end
    update_positions(system, new_positions)

end

@doc raw"""
    get_positions_cart(system::AbstractSystem, positions::Vector{<:AbstractVector{<:Real}})

    Returns positions, in cartesian coordinates, of the system's particles.

    # Returns
    - `Vector{Vector{L}}`
"""
function get_positions_cart(system::AbstractSystem)
    [a.position for a in system.particles]
end

@doc raw"""
    get_positions(system::AbstractSystem, positions::Vector{<:AbstractVector{<:Real}})

    Returns positions, in fractional coordinates, of the system's particles.

    # Returns
    - `Vector{Vector{L}}`
"""
function get_positions(system::AbstractSystem)
    pos_cart = get_positions_cart(system)
    cartesian_to_fractional(system, pos_cart)
end

@doc raw"""
    n_optimizable_coordinates(system::AbstractSystem)

Returns the number of optimizable coordinates in the system.
"""
function n_optimizable_coordinates(system::AbstractSystem)
    sum(vcat([a[:optimizable] for a in system.particles]...))
end

@doc raw"""
    set_optimizable_mask(system::AbstractSystem, mask::AbstractVector{<:AbstractVector{<:Bool}})

Sets the masks defining which coordinates of the system can be optimized. 
The mask is a list of 3-d boolean vectors, where, for each particle of the system, 
the vector specifies which of the xyz coordinates can be optimized by flagging them with `true`.
"""
function set_optimizable_mask(system::AbstractSystem, mask::AbstractVector{<:AbstractVector{<:Bool}})
    particles = [Atom(atom; optimizable=m) for (atom, m) in zip(system.particles, mask)]
    AbstractSystem(system, particles=particles)
end

@doc raw"""
    get_optimizable_mask(system::AbstractSystem) -> AbstractVector{<:AbstractVector{<:Bool}}

Returns the optimizable mask of the system (see documentation for `set_optimizable_mask`.
"""
function get_optimizable_mask(system::AbstractSystem)
    [a[:optimizable] for a in system.particles]
end

@doc raw"""
    Given a list of vectors and a mask, return the masked version in a one-dimensional list.
    This is used to concatenate the forces returned in the scfres into a single list, 
    that can be used by the optimizer.

    """
function mask_vector_list(x::AbstractVector{<:AbstractVector{<:Any}}, mask::AbstractVector{<:AbstractVector{<:Bool}})
    tmp = [element[maskelement] for (element, maskelement) in zip(x, mask)]
    # Remove empty subarrays and concatenate.
    vcat(tmp...)
end

function get_optimizable_coordinates_cart(system::AbstractSystem)
    positions = get_positions_cart(system)
    mask = get_optimizable_mask(system)
    mask_vector_list(positions, mask)
end

@doc raw"""
    Clamp given atoms if the system. Clamped atoms are fixed and their positions 
    will not be optimized. The atoms to be clamped should be given as a list of 
    indies corresponding to their positions in system.particles.

    """
function clamp_atoms(system::AbstractSystem, clamped_indexes::Union{AbstractVector{<:Integer},Nothing})
    mask = [trues(3) for _ in eachindex(system.particles)]
    if !isnothing(clamped_indexes)
        mask[clamped_indexes] = falses(3)
    end
    clamped_system = set_optimizable_mask(system, mask)
    return clamped_system
end
