# Neural network models

Using the [ASE interface](@ref) we can directly use models trained using
[SchNetPack](https://github.com/atomistic-machine-learning/schnetpack).

To use a SchNet model, load any pre-trained model into a given path you can access.
In the example below, our SchNet model is named "best_model" as is common in SchNet and
uses the the relative path.

First we load the model into an `ase` calculator and attach it to our diatomic
hydrogen molecule.
```@example spk
using PyCall

ase = pyimport("ase")
spkutils = pyimport("schnetpack.utils")
spkinterfaces = pyimport("schnetpack.interfaces")

spk_model = spkutils.load_model("../assets/schnetpack/best_model"; map_location="cpu")

h2 = ase.Atoms("H2", [(0, 0, 0), (0, 0, 0.74)])

calc = spkinterfaces.SpkCalculator(spk_model, energy="energy", forces="forces")
h2.set_calculator(calc)
```

We can obtain the energies and forces from `ase` directly in the usual way, converting
them to atomic units using [UnifulAtomic](https://github.com/sostock/UnitfulAtomic.jl).
```@repl spk
using Unitful, UnitfulAtomic;
austrip(h2.get_total_energy() * u"eV")
austrip.(h2.get_forces() .* u"eV/Å")
```

!!! warning

    Note that this is an arbitrary model not trained on H2, hence the calculation of the
    potential energy and forces most likely do not make sense.

Then, we can convert the ASE output into the format used in NonadiabaticModels,
which makes it possible to use the SchNet model e.g. for molecular dynamics calculations
within NonadiabaticModels:
```@repl spk
using NonadiabaticModels;
model = AdiabaticASEModel(h2);

r = [0 0; 0 0; 0 ustrip(auconvert(0.74u"Å"))]

potential(model, r)
derivative(model, r)
```
