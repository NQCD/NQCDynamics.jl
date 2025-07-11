```@setup logging
@info "Expanding src/NQCModels/analyticmodels.md..."
start_time = time()
```
# Analytic model library

This page plots many of the analytic models included in `NQCDynamics`.

## [`ClassicalModels`](@ref)
These models are used for classical dynamics and provide a single potential energy surface.

### [`Harmonic`](@ref)
```@example
using NQCModels, Plots
plot(-10:0.1:10, Harmonic())
```

### [`DiatomicHarmonic`](@ref)

```@example
using NQCModels, Plots

model = DiatomicHarmonic(r₀=10.0)
f(x,y) = potential(model, [x y 0])
contour(-10:0.1:10, -10:0.1:10, f, fill=true)
xlabel!("x coordinate /a₀")
ylabel!("y coordinate /a₀")
```

### [`DarlingHollowayElbow`](@ref)

```@example
using NQCModels, Plots
using NQCBase: eV_to_au

model = DarlingHollowayElbow()
V(x,z) = potential(model, [x, z])

x = range(-0.5, 3.5, length=200)
z = range(-0.5, 4.5, length=200)

plot(
    xlabel="Bond length /a₀",
    ylabel="Surface molecule distance /a₀",
    xlims=(-0.5, 3.5),
    ylims=(-0.5, 4.5)
)

contourf!(x, z, V)
```

## [`QuantumModels`](@ref QuantumModels)
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
using NQCModels, Plots

model = GatesHollowayElbow()
v1(x,z) = potential(model, [x z])[1,1]
v2(x,z) = potential(model, [x z])[2,2]

x = range(-0.5, 4.0, length=200)
z = range(-0.5, 4.0, length=200)

contour(x, z, v1, color=:blue, levels=0:0.01:0.1, label="V11", colorbar=false)
contour!(x, z, v2, color=:red, levels=0:0.01:0.1, label="V22")
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
### [`ErpenbeckThoss`](@ref)
```@example
using NQCModels, Plots

model = ErpenbeckThoss(;Γ=0.01)
ε₀(r) = potential(model, hcat(r))[1,1]
ε₁(r) = potential(model, hcat(r))[2,2]
Vₖ(r) = potential(model, hcat(r))[1,2]

r = range(1.89, 9.45, length=200)

plot(r, [ε₀, ε₁, Vₖ], label=["ε₀(r)" "ε₁(r)" "Vₖ(r)"], xlabel="r (a₀)", ylabel="V(r) (Eₕ)"; coupling=true)
```
