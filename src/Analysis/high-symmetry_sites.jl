"""
These analysis functions contain some postprocessing code to associate positions of adsorbate molecules on a surface slab with high-symmetry locations on the surface slab. 

Base definitions for fcc metal surface facets are included, but there are options to expand this further. 

**This code only works for 2D positions in X and Y, so the surface needs to lie in the XY plane.**

"""
module HighSymmetrySites


using NQCDynamics: AbstractSimulation, PeriodicCell, get_positions, NQCBase
using LinearAlgebra: norm

export SlabStructure, FCC100Sites, FCC110Sites, FCC111Sites, FCC211Sites, positions_to_category, classify_every_frame

struct SlabStructure
    adsorbate_indices::Vector{Int}
    symmetry_sites::Dict{Symbol, Vector{Vector{Float64}}}
    supercell_size::Vector{Float64}
end

const FCC100Sites = Dict(
    :top => vec([vcat(i...) for i in Iterators.product([0.0,1.0], [0.0,1.0])]),
    :bridge => [
      [0.5, 0.0],
      [0.5, 1.0],
      [0.0, 0.5],
      [1.0, 0.5],
    ],
    :fcc => [
      [0.5, 0.5],
    ],
)
const FCC110Sites = Dict(
    :top => vec([vcat(i...) for i in Iterators.product([0.0,1.0], [0.0,1.0])]),
    :long_bridge => [
      [0.5, 0.0],
      [0.5, 1.0],
    ],
    :short_bridge => [
      [0.0, 0.5],
      [1.0, 0.5],
    ],
    :center_hollow => [
      [0.5, 0.5],
    ],
    :step_hollow => [
      [0.25, 0.5],
      [0.75, 0.5],
    ],
)
const FCC111Sites = Dict(
    :top => vec([vcat(i...) for i in Iterators.product([0.0,1.0], [0.0,1.0])]),
    :fcc => [[0.33,0.33]],
    :hcp => [[0.66,0.66]],
    :bridge => [
      [0.5,0.0],
      [0.5,0.5],
      [0.0,0.5],
      [1.0,0.5],
      [0.5,1.0],
    ],
)
const FCC211Sites = Dict( # Based on Cao2018 site definitions
:step_edge => vec([vcat(i...) for i in Iterators.product([0.0,1.0], [0.0,1.0])]), # Top sites, marking desorption bridged over the step edge
:short_step => [ # This is the higher hollow site close to the step edge and would count for desorption via the step edge. Site 1 in Cao2018
  [0.125, 0.5, ] #
], 
:fcc_high => [ # Site 2 in Cao2018, equivalent to desorption not over the step edge
  [0.5, 0, ],
  [0.5, 1.0, ],
],
:fcc_low => [ # This is hollow site 3 in Cao2018, equivalent to desorption not over the step edge
  [0.75,0.0, ],
  [0.75,1.0, ],
],
:long_step => [ # This is hollow site 4, the lower one with the higher step edge, counting for desorption over the higher step.
  [0.875, 0.5, ],
], #Not really sure where b2 goes. 
)

"""
    positions_to_category(position, categories, cell::PeriodicCell; fractional::Bool = false, snap_to_site::Float64 = 0.03)

Classifies a given 2D position according to its proximity to high-symmetry surface sites.

# Arguments
- `position`: A 2-element (or longer, only first 2D used) vector representing the coordinates to classify.
- `categories`: A dictionary mapping site names (as Symbols) to lists of site coordinates.
- `cell::PeriodicCell`: The surface unit cell used to calculate distances and perform coordinate transformations.
- `fractional::Bool = false`: Whether to use fractional coordinates to determine distances to sites. If `true`, `position` is mapped directly to the sites; if `false`, converts to fractional using `cell`.
- `snap_to_site::Float64 = 0.03`: The distance tolerance (in the same units as position/cell) for snapping to a given site category.

# Returns
- The key of the closest site category (from `categories`) if within `snap_to_site` of any defined site. If not within this range or outside the cell, returns `:other`.

# Notes
- Positions are always truncated to 2D before classification.
- If the position is not inside the periodic cell, returns `:other`.
- Site assignment is based on the minimum Euclidean distance to any site in each category.
"""
function positions_to_category(
	position, 
	categories, 
	cell::PeriodicCell; 
	fractional::Bool=false, 
	snap_to_site::Float64=0.03,
)
	position = first(position,2) # Ensure we only ever carry X and Y
	category_names = collect(keys(categories))
	if !NQCBase.check_atoms_in_cell(cell, hcat(position))
		return :other
	else
		position_copy = position
		# Fractionalise positions if requested.
		fractional ? position_copy = cell.inverse * position_copy : nothing
		site_distances = Float64[]
		for site in category_names
			if fractional
				push!(site_distances, minimum([norm(position_copy - point) for point in categories[site]])) # Gets closest distance for each category
			else
				push!(site_distances, minimum([norm(position_copy - cell.vectors * point) for point in categories[site]])) # Gets closest distance for each category
			end
		end
	end
	if any(site_distances .≤ snap_to_site) # Check if any category is sufficiently close
		return category_names[argmin(site_distances)] # Return the classified category
	else
		return :other
	end
end

"""
    classify_every_frame(
        trajectory,
        cell::PeriodicCell,
        slab::SlabStructure,
        fractional::Bool = false,
        snap_to_site::Float64 = 0.03,
    )

Classifies the adsorption site of adsorbate atoms for every frame in a trajectory.

# Arguments
- `trajectory`: Vector{DynamicsVariables type} of structures to classify. 
- `cell::PeriodicCell`: The periodic cell representation of the simulation.
- `slab::SlabStructure`: Slab structure object containing adsorbate atom indices, supercell size, and symmetry site information.
- `fractional::Bool = false`: If `true`, uses fractional coordinates for site matching; otherwise, uses Cartesian coordinates.
- `snap_to_site::Float64 = 0.03`: Distance tolerance for assigning a position to a high-symmetry site.

# Returns
- A vector of vectors, where each inner vector contains the classified site (as a `Symbol`) for an adsorbate atom at one frame.

# Notes
- Each adsorbate atom is assigned to a symmetry site or `:other` for each trajectory frame.
- Positions are wrapped back into the primitive cell before site classification.
- Uses `positions_to_category` to perform the categorization for each atom in each frame.
"""
function classify_every_frame(
	trajectory, 
    cell::PeriodicCell,
    slab::SlabStructure;
	fractional::Bool = false,
	snap_to_site::Float64 = 0.03,
)
	# Determine a 2D representation of the unit cell first
    fullsize_cell = cell.vectors
	primitive_cell = PeriodicCell(fullsize_cell ./ slab.supercell_size)
	primitive_cell_2d = PeriodicCell(getindex(fullsize_cell ./ slab.supercell_size, 1:2, 1:2))

	site_associations = [Symbol[] for _ in eachindex(slab.adsorbate_indices)] # Generate site association vector containing info for each adsorbate atom. 
	@inbounds for trajectory_frame in eachindex(trajectory)
        @inbounds for (output_list_index, atom_index) in enumerate(slab.adsorbate_indices)
            position_H = hcat(
                get_positions(
                    trajectory[trajectory_frame]
                )[:, atom_index]
            )
            # Move into primitive unit cell
            NQCBase.apply_cell_boundaries!(primitive_cell, position_H)
            position_H = vec(position_H)
            # Then categorise
            category = positions_to_category(
                position_H,
                slab.symmetry_sites,
                primitive_cell_2d,
                fractional=fractional,
                snap_to_site=snap_to_site
            )
            push!(site_associations[output_list_index], category)
        end
	end
	return site_associations
end


end