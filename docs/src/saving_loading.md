```@setup logging
@info "Expanding src/saving_loading.md..."
start_time = time()
```
# [Saving and loading](@id saving-and-loading)

If you would like to split your workflow into multiple scripts
(e.g. separately generating initial conditions and running dynamics)
it is necessary to be able to store intermediate data in files.

When this data is atomic configurations or trajectories, it can be useful to use
standard file formats such as those mentioned in the [Atoms section](@ref atoms) previously.
However, often it is more convenient to directly save and load Julia objects between sessions.
For this purpose, we recommend using [FileIO.jl](https://github.com/JuliaIO/FileIO.jl) with [JLD2.jl](https://github.com/JuliaIO/JLD2.jl).

!!! note

    JLD2 can be used independently of FileIO. However, FileIO provides a unified interface
    for many file types and allows you to save data to lots of formats with consistent
    syntax.

As a simple example, suppose that we want the same system parameters across multiple scripts:

```@example saving
using NQCDynamics

atoms = Atoms([:H, :H, :C, :C])
cell = PeriodicCell([10.0 0 0; 0 10.0 0; 0 0 10.0])
model = Harmonic()
nothing # hide
```

Instead of redefining these in every script, we can save them to a file, then load them back in whenever we need them
using FileIO.

This creates a file `"parameters.jld2"` containing all of our system parameters:
```@example saving
using FileIO
save("parameters.jld2", Dict("atoms"=>atoms, "cell"=>cell, "model"=>model))
```

In a separate Julia session we can re-load these parameters.
As detailed in the [JLD2 documentation](https://github.com/JuliaIO/JLD2.jl) we
can select the data to load by specifying extra arguments to `load`.
```@repl saving2
using NQCDynamics, FileIO

parameters = load("parameters.jld2")
atoms = load("parameters.jld2", "atoms")
cell = load("parameters.jld2", "cell")
model = load("parameters.jld2", "model")
```

[JLD2](https://github.com/JuliaIO/JLD2.jl) is compatible with any Julia type so it widely
usable for most of the types you encounter is NQCDynamics.jl and across all Julia packages.
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
