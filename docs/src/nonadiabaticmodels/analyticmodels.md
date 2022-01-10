# Analytic model library

This page plots many of the analytic models included in `NQCDynamics`.

!!! tip 
    To produce the plots we use two of Julia's plotting options [Plots](http://docs.juliaplots.org/latest/)
    and [Makie](https://makie.juliaplots.org/dev/).
    Plots has a mature recipe system that allows us to define custom plots for the 1D models
    but we use Makie to produce the more complex images.
    Each of these has their pros and cons and if you are interested in producing plots
    using Julia you should visit their documentation to decide which is best for you.

## [`AdiabaticModels`](@ref NQCModels.AdiabaticModels)
These models are used for classical dynamics and provide a single potential energy surface.

### [`Harmonic`](@ref)
```@example
using NQCModels, Plots
plot(-10:0.1:10, Harmonic())
```

### [`DiatomicHarmonic`](@ref)

```@example
using NQCModels
using Plots

model = DiatomicHarmonic(r₀=10.0)
f(x,y) = potential(model, [x y 0])
contour(-10:0.1:10, -10:0.1:10, f, fill=true)
xlabel!("x coordinate /a₀")
ylabel!("y coordinate /a₀")
```

### [`DarlingHollowayElbow`](@ref)

```@example
using NQCModels
using NQCBase: eV_to_au
using CairoMakie

model = DarlingHollowayElbow()
V(x,z) = potential(model, [x, z])

x = range(-0.5, 3.5, length=200)
z = range(-0.5, 4.5, length=200)

f = Figure()
ax = Axis(f[1,1], xlabel="Bond length /a₀", ylabel="Surface molecule distance /a₀")

levels = eV_to_au.([-0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1])
contourf!(ax, x, z, V, levels=levels)
contour!(ax, x, z, V, levels=levels, color=:black)
Colorbar(f[1,2], limits=(-0.1, 1.1))
xlims!(-0.5, 3.5)
ylims!(-0.5, 4.5)
f
```

## [`DiabaticModels`](@ref NQCModels.DiabaticModels)
These models define a Hermitian potential operator in a diabatic basis.
These can be used for various forms of nonadiabatic dynamics.

### [`TullyModelOne`](@ref)
```@example
using NQCModels, Plots
plot(-10:0.1:10, TullyModelOne(); coupling=true)
```
### [`TullyModelTwo`](@ref)
```@example
using NQCModels, Plots
plot(-10:0.1:10, TullyModelTwo(); coupling=true)
```
### [`TullyModelThree`](@ref)
```@example
using NQCModels, Plots
plot(-10:0.1:10, TullyModelThree(); coupling=true)
```
### [`ThreeStateMorse`](@ref)
```@example
using NQCModels, Plots
plot(2:0.01:12, ThreeStateMorse(), ylims=(0, 0.06), coupling=true)
```
### [`OuyangModelOne`](@ref)
```@example
using NQCModels, Plots
plot(-10:0.1:10, OuyangModelOne())
```
### [`DoubleWell`](@ref)
```@example
using NQCModels, Plots
plot(-5:0.1:5, DoubleWell(); coupling=true)
```
### [`GatesHollowayElbow`](@ref)
```@example
using NQCModels
using CairoMakie

model = GatesHollowayElbow()
v1(x,z) = potential(model, [x z])[1,1]
v2(x,z) = potential(model, [x z])[2,2]
coupling(x,z) = potential(model, [x z])[1,2]

x = range(-0.5, 4.0, length=200)
z = range(-0.5, 4.0, length=200)

f = Figure()
ax = Axis(f[1,1], xlabel="x coordinate", ylabel="z coordinate")

contour!(ax, x, z, coupling, color=:black, levels=10, label="V12")
contour!(ax, x, z, v1, color=:blue, levels=0:0.01:0.1, label="V11")
contour!(ax, x, z, v2, color=:red, levels=0:0.01:0.1, label="V22")
axislegend(ax)
xlims!(-0.5, 4.0)
ylims!(-0.5, 4.0)
f
```
