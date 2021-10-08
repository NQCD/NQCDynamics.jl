# ASE interface

The easiest way to obtain potentials and forces from established codes is to directly
use the interfaces implemented in [ASE](https://wiki.fysik.dtu.dk/ase/).

We provide the [`AdiabaticASEModel`](@ref) which wraps an ASE atoms object and its
associated calculator to implement the required [`potential`](@ref) and
[`derivative`](@ref) functions.

!!! note

    The interface works by directly calling the relevant Python functions using
    [PyCall](https://github.com/JuliaPy/PyCall.jl).

!!! note

    To use PyCall, you must make sure that the python version contains all relevant packages, such as ase. Alternatively, you can add the Pkg "conda", which installs miniconda private to julia. With this, you can then install necessary packages, such as ase.
    
## Example

First, it is necessary to import `ase` and create the `ase.Atoms` object and attach
the desired calculator. This works exactly as in Python:
```@example ase
using PyCall

ase = pyimport("ase")
emt = pyimport("ase.calculators.emt")

h2 = ase.Atoms("H2", [(0, 0, 0), (0, 0, 0.74)])
h2.calc = emt.EMT()
nothing # hide
```

Next, the [`AdiabaticASEModel`](@ref) is created by passing the `ase.Atoms` object directly
to the model:
```@repl ase
using NonadiabaticModels
model = AdiabaticASEModel(h2)
```
Now the model can be used in the same way as any of the previously introduced
analytic models.
```@repl ase
potential(model, rand(3, 2))
derivative(model, rand(3, 2))
```

!!! tip 

    In theory, this should work with any of the ASE calculators that correctly implement
    the `get_potential_energy` and `get_forces` functions.
