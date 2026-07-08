# [Simulations outputs and analysis functions](@id output_and_analysis)

NQCDynamics.jl is able to output a variety of common observables when running simulations with the `run_dynamics` function. Further output functions can be implemented easily as well. 

In this section you will find an overview of all available output types, as well as explanations of some common analysis functions in the realm of surface chemistry which we have implemented in the package. 

## DynamicsOutputs

In many examples within this documentation, you have seen calls to `run_dynamics`:

```julia
ensemble = run_dynamics(sim, tspan, distribution;
    dt=0.1u"fs",
    output=OutputPosition,
    trajectories=20
)
```

Within `run_dynamics`, the `output` argument specifies the desired output values for each trajectory. `output` can either be given as a single function, or as a tuple of multiple output functions, for example:

```julia
output=OutputPosition # or
output=(OutputPosition, OutputVelocity, OutputKineticEnergy)
```

The "ensemble" variable that `run_dynamics` outputs to is a vector of dictionaries, where each dictionary in the vector relates to one of the trajectories. By indexing the `ensemble` vector all outputs for that trajectory will be returned.

```
ensemble[3][:OutputPosition] # will output the positions at all timesteps in trajectory 3
```

Every output type is a function which can use the [`DynamicsVariables`](@ref) and [`Simulation`](@ref) values of the respective trajectory, allowing you to create custom output types of any kind. See the [developer documentation] for more information on how to implement a custom output type. 

You can find an overview of all available output types in the [`DynamicsOutputs`](@ref NQCDynamics.DynamicsOutputs) API. 

## Analysis functions

The [Analysis](@ref Analysis) submodule in NQCDynamics contains functions commonly used in the analysis of trajectories to make the analysis of existing trajectories easier. 
Ideally, most observable quantities could be implemented with a combination of [`DynamicsOutputs`](@ref NQCDynamics.DynamicsOutputs) and `Reduction` types, however we might want to data from existing ensemble simulations where re-running the entire set of trajectories is impractical. 

As a result, most functions in the `Analysis` submodule are also implemented as a `DynamicsOutput`. 

### Convenient functions for periodic structures
[`NQCDynamics.Structure`](@ref Structure) contains several useful functions for periodic structures, such as `pbc_distance, pbc_center_of_mass`. 

These functions take into account periodic copies of the atoms in question, returning the respective values for the closest set of periodic copies. 

### Analysis of diatomic molecules
[`NQCDynamics.Analysis.Diatomic`](@ref Analysis)

### Surface site classification
[`NQCDynamics.Analysis.HighSymmetrySites`](@ref Analysis.HighSymmetrySites) provides functionality for classifying adsorbate positions on periodic surface slabs according to high-symmetry surface sites.

This module is particularly useful for analyzing trajectories of molecules on metal surfaces, where you want to track whether adsorbates are located at top, bridge, hollow, or other characteristic sites.

**Important:** This code only works for 2D positions in X and Y, so the surface must lie in the XY plane.

#### Predefined surface facets

The module includes predefined site definitions for common FCC metal surface facets:
- `FCC100Sites` - (100) surface with top, bridge, and hollow sites
- `FCC110Sites` - (110) surface with top, long bridge, short bridge, center hollow, and step hollow sites
- `FCC111Sites` - (111) surface with top, bridge, and hollow sites
- `FCC211Sites` - (211) stepped surface with step edge, short step, fcc high, fcc low, and long step sites

#### Basic usage

To classify adsorbate positions, you need to:
1. Define a `SlabStructure` containing:
   - `adsorbate_indices`: indices of atoms to track
   - `symmetry_sites`: dictionary mapping site names to their fractional coordinates
   - `supercell_size`: size of the simulation cell relative to the primitive cell
2. Use `positions_to_category` to classify a single position, or
3. Use `classify_every_frame` to analyze an entire trajectory

#### Example

```julia
using NQCDynamics

# Define the slab structure for an FCC(100) surface
cell = PeriodicCell(diagm([10.0, 10.0, 20.0]))  # 10×10 Å surface cell
slab = Analysis.HighSymmetrySites.SlabStructure(
    [1, 2],  # Track atoms 1 and 2
    Analysis.HighSymmetrySites.FCC100Sites,  # Use FCC(100) site definitions
    [2.0, 2.0, 1.0]  # 2×2 supercell in XY
)

# Classify positions from a trajectory
# trajectory should be a vector of DynamicsVariables
sites = Analysis.HighSymmetrySites.classify_every_frame(
    trajectory,
    cell,
    slab,
    snap_to_site=0.3  # Distance tolerance in Angstroms
)

# sites[i] contains the classified site for adsorbate i at each frame
# Each site is a Symbol like :top, :bridge, :fcc, or :other
```

For more details, see the [API documentation](@ref Analysis.HighSymmetrySites).
