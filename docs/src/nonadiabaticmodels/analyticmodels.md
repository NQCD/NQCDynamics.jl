# Analytic model library

This page plots many of the analytic models included in `NonadiabaticMolecularDynamics`.

!!! tip 
    To produce the plots we use two of Julia's plotting options [`Plots`](http://docs.juliaplots.org/latest/)
    and [`Makie`](https://makie.juliaplots.org/dev/).
    Plots has a mature recipe system that allows us to define custom plots for the 1D models
    but we use Makie to produce the more complex images.

```@setup model
using NonadiabaticModels
using Plots: plot, contour
```

## [`AdiabaticModels`](@ref NonadiabaticModels.AdiabaticModels)
These models are used for classical dynamics and provide a single potential energy surface.

### [`DiatomicHarmonic`](@ref)

```@example model
model = DiatomicHarmonic(r₀=1.0)
f(x,y) = potential(model, [x y])[1]
contour(-10:0.1:10, -10:0.1:10, f, fill=true)
```

### [`DarlingHollowayElbow`](@ref)

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

## [`DiabaticModels`](@ref NonadiabaticModels.DiabaticModels)
These models define a Hermitian potential operator in a diabatic basis.
These can be used for various forms of nonadiabatic dynamics.

### [`TullyModelOne`](@ref)
```@example model
plot(-10:0.1:10, TullyModelOne())
```
### [`TullyModelTwo`](@ref)
```@example model
plot(-10:0.1:10, TullyModelTwo())
```
### [`TullyModelThree`](@ref)
```@example model
plot(-10:0.1:10, TullyModelThree())
```
### [`ThreeStateMorse`](@ref)
```@example model
plot(2:0.01:12, ThreeStateMorse(), ylims=(0, 0.06))
```
### [`OuyangModelOne`](@ref)
```@example model
plot(-10:0.1:10, OuyangModelOne())
```
### [`DoubleWell`](@ref)
```@example model
plot(-5:0.1:5, DoubleWell())
```
### [`GatesHollowayElbow`](@ref)
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
