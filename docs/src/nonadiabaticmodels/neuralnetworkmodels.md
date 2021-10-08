# Neural network models

Using the [ASE interface](@ref) we can directly use models trained using
[SchNetPack](https://github.com/atomistic-machine-learning/schnetpack).

To use a SchNet model, please load any pre-trained model into a given path you can access.
Here, our SchNet model is named "best_model" as is common in SchNet and
provide the relative path.

First we can load the model into an `ase` calculator and attach it to our diatomic
hydrogen molecule.
```@example spk
using PyCall

ase = pyimport("ase")
spk = pyimport("schnetpack")

spk_model = spk.utils.load_model("../assets/schnetpack/best_model")

h2 = ase.Atoms("H2", [(0, 0, 0), (0, 0, 0.74)])

calc = spk.interfaces.SpkCalculator(model, energy="energy", forces="forces")
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

Then we can obtain the same numbers using the NonadiabaticModels interface:
```@repl spk
using NonadiabaticModels;

model = AdiabaticASEModel(spk_model)

r = [0 0; 0 0; 0 auconvert(0.74u"Å")]

potential(model, r)
derivative(model, r)
```
