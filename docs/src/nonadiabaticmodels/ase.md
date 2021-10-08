# ASE interface

The easiest way to obtain potentials and forces from established codes is to directly
use the interfaces implemented in [ASE](https://wiki.fysik.dtu.dk/ase/).

We provide the [`AdiabaticASEModel`](@ref) which wraps an ASE atoms object and its
associated calculator to implement the required [`potential`](@ref) and
[`derivative`](@ref) functions.

!!! note

    The interface works by directly calling the relevant Python functions using
    [PyCall](https://github.com/JuliaPy/PyCall.jl).

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
    the `get_potential_energy` and `get_forces` functions. For isnstance, you can install SchNetPack (SPK) in your python version and import it via "spk=pyimport("schnetpack"). You then have to load a model and use this as the calculator. 


To use a SchNet model, please load any pre-trained model into a given path you can access. Here, we assume that the SchNet model is named "best_model" as usually in SchNet and is in the same path where we are executing the code.

```@repl spk
using PyCall

ase = pyimport("ase")
spk = pyimport("schnetpack")

#load model
model = spk.utils.load_model("best_model")

#define ase-atoms type
h2 = ase.Atoms("H2", [(0, 0, 0), (0, 0, 0.74)])

#set SPK calculator
calc = spk.interfaces.SpkCalculator(model, energy="energy", forces="forces")
h2.set_calculator(calc)

#predict energy and forces

h2.get_total_energy()
h2.get_forces()```

!!! note
Note that this is an arbitrary model not trained on H2, hence the calculation of the potential energy and forces most likely do not make sense.
