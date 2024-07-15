```@setup logging
@info "Expanding src/NQCModels/neuralnetworkmodels.md..."
start_time = time()
```

# MLIP Examples

## SchNet (SchNetPack) models

Using the [ASE interface](@ref) we can directly use models trained using
[SchNetPack](https://github.com/atomistic-machine-learning/schnetpack).

!!! warning

    The examples on this page do not run during the documentation build due to `schnetpack`
    causing segfaults when installed in the build environment.
    The causes of this is not currently clear but we have temporarily disabled these examples
    in the build.

    However, the examples should still be correct and you are welcome to try them with
    your own schnetpack trained models.

To use a SchNet model, please load any pre-trained model into a given path you can access.
Here, our SchNet model is named "best_model" as is common in SchNet and
provide the relative path.

First we load the model into an `ase` calculator and attach it to our diatomic
hydrogen molecule.

```julia
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
them to atomic units using [UnitfulAtomic](https://github.com/sostock/UnitfulAtomic.jl).

```julia-repl
using Unitful, UnitfulAtomic;
austrip(h2.get_total_energy() * u"eV")
austrip.(h2.get_forces() .* u"eV/Å")
```

!!! warning

    Note that this is an arbitrary model not trained on H2, hence the calculation of the
    potential energy and forces most likely do not make sense.

Then, we can convert the ASE output into the format used in NQCModels,
which makes it possible to use the SchNet model e.g. for molecular dynamics calculations
within NQCDynamics.jl:

```julia-repl
using NQCModels;
model = AdiabaticASEModel(h2);

r = [0 0; 0 0; 0 ustrip(auconvert(0.74u"Å"))]

potential(model, r)
derivative(model, r)
```

## MACE models

!!! warning

      A more direct interface for the use of MACE models is now available. @Alex link to example page and API docs here.

      The following example is still valid but the new interface is recommended for new projects.

The following example shows how to connect a trained [MACE](https://github.com/ACEsuit/mace) model to NQCDynamics.jl using the [ASE interface](@ref).

```julia
using PyCall, NQCModels

ase_io = pyimport("ase.io") # make sure ase is installed in your environment first
mace_module = pyimport("mace.calculators") # make sure mace is installed in your environment first

structure = ase_io.read("starting_structure.xyz")
mace_calculator = mace_module.MACECalculator(
    "path/to/model.model",
    device = "cpu", # or cuda to run on GPU
    default_dtype = "float32" # if this was set during training as well
)

structure.set_calculator(mace_calculator)
model = AdiabaticASEModel(structure)
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
