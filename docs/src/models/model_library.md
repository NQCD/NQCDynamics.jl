# Analytic model library

Provided by the package are a few analytic models that can be used to test new
algorithms. Here you can find descriptions of the models and examples of their usage. 

!!! tip

    In some of the examples we use `Symbolics.jl` to show the functional form of the output,
    this can be useful for checking the model has been defined correctly.
    The models have been defined in a type-agnostic way such that they should be easily
    extensible with custom types.

```@setup model
using NonadiabaticModels
using Plots: plot, contour
```

## Adiabatic models
These models are used for classical dynamics and provide a single potential energy surface.
```@docs
NonadiabaticModels.Free
NonadiabaticModels.Harmonic
NonadiabaticModels.JuLIPModel
```

```@docs
NonadiabaticModels.DiatomicHarmonic
```
```@example model
model = DiatomicHarmonic(r₀=1.0)
f(x,y) = potential(model, [x y])[1]
contour(-10:0.1:10, -10:0.1:10, f, fill=true)
```

```@docs
NonadiabaticModels.DarlingHollowayElbow
```
```@example
using NonadiabaticModels # hide
using NonadiabaticDynamicsBase: eV_to_au
using CairoMakie

model = DarlingHollowayElbow()
V(x,z) = potential(model, [x z])[1]

x = range(-0.5, 3.5, length=200)
z = range(-0.5, 4.5, length=200)

fig = Figure()

levels = eV_to_au.([-0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1])
contourf(fig[1,1], x, z, V, levels=levels)
contour!(fig[1,1], x, z, V, levels=levels, color=:black)
Colorbar(fig[1,2], limits=(-0.1, 1.1))
xlims!(-0.5, 3.5)
ylims!(-0.5, 4.5)
fig
```

```@docs
NonadiabaticModels.DebyeBosonBath
```
## Diabatic models
These models define a Hermitian potential operator in a diabatic basis.
These can be used for various forms of nonadiabatic dynamics.
```@docs
NonadiabaticModels.TullyModelOne
```
```@example model
plot(-10:0.1:10, TullyModelOne())
```
```@docs
NonadiabaticModels.TullyModelTwo
```
```@example model
plot(-10:0.1:10, TullyModelTwo())
```
```@docs
NonadiabaticModels.ThreeStateMorse
```
```@example model
plot(2:0.01:12, ThreeStateMorse(), ylims=(0, 0.06))
```
```@docs
NonadiabaticModels.OuyangModelOne
```
```@example model
plot(-10:0.1:10, OuyangModelOne())
```
```@docs
NonadiabaticModels.DoubleWell
```
```@example model
plot(-5:0.1:5, DoubleWell())
```
```@docs
NonadiabaticModels.GatesHollowayElbow
```
```@example
using NonadiabaticModels # hide
using CairoMakie

model = GatesHollowayElbow()
v1(x,z) = potential(model, [x z])[1,1]
v2(x,z) = potential(model, [x z])[2,2]
coupling(x,z) = potential(model, [x z])[1,2]

x = range(-0.5, 4.0, length=200)
z = range(-0.5, 4.0, length=200)

fig = Figure()

contour(fig[1,1], x, z, coupling, color=:black, levels=10)
contour!(fig[1,1], x, z, v1, color=:blue, levels=0:0.01:0.1)
contour!(fig[1,1], x, z, v2, color=:red, levels=0:0.01:0.1)
xlims!(-0.5, 4.0)
ylims!(-0.5, 4.0)
fig
```

```@docs
NonadiabaticModels.DebyeSpinBoson
```