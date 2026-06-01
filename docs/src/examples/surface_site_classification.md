```@setup logging
@info "Expanding src/examples/surface_site_classification.md..."
start_time = time()
```
# Surface site classification for adsorbates

This example demonstrates how to use the `HighSymmetrySites` analysis module to classify the positions of adsorbate molecules on a metal surface according to high-symmetry surface sites.

The `HighSymmetrySites` module is useful for analyzing molecular dynamics trajectories on periodic surfaces, particularly for studying adsorption, desorption, and surface diffusion processes.

## Setting up the system

We'll create a simple example with a hydrogen molecule on an FCC(100) metal surface.
First, let's set up the basic simulation cell and atoms.

```@example surface_sites
using NQCDynamics
using LinearAlgebra

# Create a periodic cell for the surface (10×10 Å in XY, 20 Å in Z)
cell = PeriodicCell(diagm([10.0, 10.0, 20.0]))

# Create atoms: 4 surface atoms and 2 hydrogen atoms (adsorbate)
atoms = Atoms([:Au, :Au, :Au, :Au, :H, :H])
nothing # hide
```

## Defining the slab structure

To use the site classification functionality, we need to define a `SlabStructure` that contains:
- The indices of adsorbate atoms to track
- The symmetry site definitions for the surface
- The supercell size relative to the primitive unit cell

For an FCC(100) surface, predefined site definitions are available via `FCC100Sites`.

```@example surface_sites
using NQCDynamics.Analysis.HighSymmetrySites

# Define which atoms are adsorbates (indices 5 and 6 for the H atoms)
adsorbate_indices = [5, 6]

# The simulation cell is a 2×2 supercell of the primitive FCC(100) unit cell
supercell_size = [2.0, 2.0, 1.0]

# Create the slab structure
slab = SlabStructure(
    adsorbate_indices,
    FCC100Sites,  # Predefined sites: :top, :bridge, :fcc
    supercell_size
)
nothing # hide
```

## Classifying a single position

Let's classify a single position to see which surface site it corresponds to.
The `positions_to_category` function determines the closest high-symmetry site.

```@example surface_sites
# Position at the origin (a top site in FCC100)
position_top = [0.0, 0.0, 5.0]  # Only X,Y coordinates are used
site_top = positions_to_category(
    position_top,
    FCC100Sites,
    PeriodicCell(diagm([5.0, 5.0]))  # Primitive cell (half of our supercell)
)
println("Position [0.0, 0.0] is classified as: ", site_top)

# Position at the bridge site
position_bridge = [2.5, 0.0, 5.0]  # At the bridge between two top sites
site_bridge = positions_to_category(
    position_bridge,
    FCC100Sites,
    PeriodicCell(diagm([5.0, 5.0])),
    snap_to_site=0.5  # Distance tolerance in Angstroms
)
println("Position [2.5, 0.0] is classified as: ", site_bridge)

# Position at the hollow site
position_hollow = [2.5, 2.5, 5.0]  # At the hollow (fcc) site
site_hollow = positions_to_category(
    position_hollow,
    FCC100Sites,
    PeriodicCell(diagm([5.0, 5.0])),
    snap_to_site=0.5
)
println("Position [2.5, 2.5] is classified as: ", site_hollow)
```

## Analyzing a trajectory

To analyze an entire trajectory, we can use `classify_every_frame`.
This function processes all frames and returns the site classification for each adsorbate at each timestep.

```@example surface_sites
# Create a simple model for demonstration
model = Free(3)  # Free particle model
sim = Simulation(atoms, model; cell=cell)

# Create some example trajectory data with different positions
positions = [
    [0.0 0.0 0.0 0.0 0.1 0.1;    # Frame 1: H atoms near top site
     0.0 0.0 0.0 0.0 0.1 0.1;
     10.0 10.0 10.0 10.0 10.0 10.0],
    [0.0 0.0 0.0 0.0 2.5 2.5;    # Frame 2: H atoms near bridge site
     0.0 0.0 0.0 0.0 0.1 0.1;
     10.0 10.0 10.0 10.0 10.0 10.0],
    [0.0 0.0 0.0 0.0 2.4 2.4;    # Frame 3: H atoms near hollow site
     0.0 0.0 0.0 0.0 2.6 2.6;
     10.0 10.0 10.0 10.0 10.0 10.0],
]

# Create DynamicsVariables for each frame
trajectory = [DynamicsVariables(sim, rand(3,6), pos) for pos in positions]

# Classify all frames
site_classifications = classify_every_frame(
    trajectory,
    cell,
    slab,
    snap_to_site=0.5  # Tolerance for site assignment
)

# Display results
println("\nTrajectory analysis:")
for (atom_idx, atom_sites) in enumerate(site_classifications)
    println("Adsorbate atom $(adsorbate_indices[atom_idx]):")
    for (frame, site) in enumerate(atom_sites)
        println("  Frame $frame: $site")
    end
end
```

## Available surface facets

The `HighSymmetrySites` module includes predefined site definitions for several common FCC metal surface facets:

### FCC(100) - `FCC100Sites`
- `:top` - Atop sites (directly above a surface atom)
- `:bridge` - Bridge sites (between two surface atoms)
- `:fcc` - Four-fold hollow site

### FCC(110) - `FCC110Sites`
- `:top` - Atop sites
- `:long_bridge` - Long bridge sites
- `:short_bridge` - Short bridge sites
- `:center_hollow` - Center hollow sites
- `:step_hollow` - Step hollow sites

### FCC(111) - `FCC111Sites`
- `:top` - Atop sites
- `:bridge` - Bridge sites
- `:hollow` - Three-fold hollow sites

### FCC(211) - `FCC211Sites`
- `:step_edge` - Step edge sites
- `:short_step` - Short step sites
- `:fcc_high` - High FCC hollow sites
- `:fcc_low` - Low FCC hollow sites
- `:long_step` - Long step sites

## Custom site definitions

You can also define custom site positions for other surface structures:

```@example surface_sites
# Define custom sites (in fractional coordinates)
custom_sites = Dict(
    :custom_site_A => [[0.25, 0.25], [0.75, 0.75]],
    :custom_site_B => [[0.25, 0.75], [0.75, 0.25]],
)

# Use with your own SlabStructure
custom_slab = SlabStructure(
    [5, 6],  # Adsorbate indices
    custom_sites,
    [2.0, 2.0, 1.0]
)
nothing # hide
```

## Tips for usage

1. **Distance tolerance**: The `snap_to_site` parameter controls how close a position must be to a defined site to be classified. Positions further away are classified as `:other`.

2. **Fractional coordinates**: Set `fractional=true` if your positions are already in fractional coordinates relative to the unit cell.

3. **Periodic boundaries**: The module automatically handles periodic boundary conditions, so adsorbates near cell edges are correctly classified.

4. **2D only**: This module only considers X and Y coordinates. The surface must lie in the XY plane.

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
