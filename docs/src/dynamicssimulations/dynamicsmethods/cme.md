```@setup logging
@info "Expanding src/dynamicssimulations/dynamicsmethods/cme.md..."
start_time = time()
```
# [Classical Master Equation (CME/BCME)](@id cme-dynamics)

The Classical Master Equation (CME) and its Broadened extension (BCME) are surface hopping methods designed specifically for two-state open quantum systems, particularly the Newns-Anderson (Anderson-Holstein) model [Dou2015](@cite)[Dou2020](@cite).
These methods describe the dynamics of an adsorbate interacting with a metallic substrate, where the adsorbate electronic state can be either empty (ground state) or occupied (excited state), and transitions between these states are driven by coupling to the metal electrons.

## CME

In CME, the classical nuclear coordinates evolve on either the neutral (empty) or charged (occupied) diabatic potential energy surface.
The hopping probability between the two states is derived from the Fermi golden rule and depends on the coupling strength ``\Gamma``, the energy gap ``\Delta V`` between the two states, and the Fermi-Dirac distribution function ``f``:

```math
P_{1 \to 2} = \Gamma f(\Delta V) \, dt
```
```math
P_{2 \to 1} = \Gamma (1 - f(\Delta V)) \, dt
```

where ``\Gamma = 2\pi V_{12}^2`` is the state-coupling-derived broadening, and ``\Delta V = V_{22} - V_{11}`` is the energy gap between the charged and neutral diabatic surfaces.
These transition probabilities ensure detailed balance, correctly describing charge transfer processes in the wide-band limit.

Velocity rescaling is not performed in CME; when a hop occurs the trajectory simply switches potential energy surface without adjusting the nuclear velocities.

## BCME

BCME (Broadened CME) extends CME by incorporating quantum broadening effects into the force calculation [Dou2016](@cite)[Dou2020](@cite).
Rather than simply using the Fermi-Dirac distribution evaluated at the current energy gap, BCME includes the effect of the finite electronic lifetime of the adsorbate level through a convolution of the Fermi function with a Lorentzian spectral density:

```math
n(\Delta V) = \int A(\epsilon) f(\epsilon) \, d\epsilon
```

where ``A(\epsilon) = \frac{1}{\pi} \frac{\Gamma/2}{(\epsilon - \Delta V)^2 + (\Gamma/2)^2}`` is the Lorentzian spectral function.

This broadening modifies both the population and the forces acting on the nuclei, giving a more accurate description of the system dynamics in the strong coupling or high temperature regime.

!!! note

    Both CME and BCME use an adaptive time stepping integrator by default.
    If you encounter numerical instabilities or unexpected results, you can increase the stability by setting tighter tolerances with the `abstol` and `reltol` keyword arguments in `run_dynamics`, or disable adaptive stepping with `adaptive=false`.

## Example

Here we demonstrate a simulation using CME for the Anderson-Holstein model.
The parameters are taken from [Dou2015](@cite), where the model describes a harmonic oscillator coupled to a two-level system representing an adsorbate level on a metal surface.

```@example cme
using NQCDynamics
using Distributions: Normal

ħω = 0.003
Γ = 0.003
kT = 0.03
g = 0.005
Ed = 0.0
atoms = Atoms(1/ħω)

struct AndersonHolsteinModel{T} <: NQCModels.QuantumModels.QuantumModel
    ħω::T
    Ed::T
    g::T
    Γ::T
end

NQCModels.nstates(::AndersonHolsteinModel) = 2
NQCModels.ndofs(::AndersonHolsteinModel) = 1

function NQCModels.potential!(model::AndersonHolsteinModel, V::Hermitian, r::AbstractMatrix)
    (;ħω, Ed, g) = model
    potential = ħω*first(r)^2/2
    V.data .= [potential sqrt(model.Γ/2π); sqrt(model.Γ/2π) potential + Ed + sqrt(2)*g*first(r)]
    return nothing
end

function NQCModels.derivative!(model::AndersonHolsteinModel, D::AbstractMatrix{<:Hermitian}, r::AbstractMatrix)
    (;ħω, g) = model
    D[1].data .= [ħω*first(r) 0; 0 ħω*first(r) + sqrt(2)*g]
    return nothing
end

model = AndersonHolsteinModel(ħω, Ed, g, Γ)
```

A `Simulation` is set up by specifying the method type `CME` along with a temperature parameter:

```@example cme
sim = Simulation{CME}(atoms, model; temperature=kT)
```

The initial distribution samples velocities from a Boltzmann distribution and positions from the appropriate thermal distribution of the harmonic oscillator.
The electronic state is initialised in the neutral (ground) state using `PureState(1, Diabatic())`:

```@example cme
βω = ħω / kT
σ = sqrt(1 / βω)
r_dist = hcat(Normal(0.0, σ))
v_dist = VelocityBoltzmann(kT, atoms.masses, (1,1))
distribution = DynamicalDistribution(v_dist, r_dist, (1,1)) * PureState(1, Diabatic())
```

We can then run an ensemble of trajectories and track the kinetic energy and electronic state over time:

```@example cme
using Statistics: mean

output = run_dynamics(sim, (0.0, 200/Γ), distribution;
    trajectories=100,
    output=(OutputKineticEnergy, OutputDiscreteState),
    abstol=1e-8, reltol=1e-8,
    saveat=10/Γ
)
```

The average kinetic energy over the ensemble, normalised by ``kT``, shows the phonon relaxation dynamics:

```@example cme
using Plots

times = output[1][:Time] .* Γ
avg_ke = mean(o[:OutputKineticEnergy] for o in output) ./ kT
plot(times, avg_ke; xlabel="t·Γ", ylabel="⟨Ekin⟩ / kT", label="CME")
```

## BCME Example

For BCME, the setup is very similar, but the `bandwidth` parameter controlling the integration range of the Lorentzian convolution must be specified.
The bandwidth should be large enough to capture the full spectral weight of the Lorentzian:

```@example bcme
using NQCDynamics
using Distributions: Normal

ħω = 0.003
Γ = 0.1
kT = 0.01
g = 0.015
Ed = g^2 / ħω
atoms = Atoms(1/ħω)
x1 = -sqrt(2) * g / ħω

struct AndersonHolsteinModel{T} <: NQCModels.QuantumModels.QuantumModel
    ħω::T
    Ed::T
    g::T
    Γ::T
end

NQCModels.nstates(::AndersonHolsteinModel) = 2
NQCModels.ndofs(::AndersonHolsteinModel) = 1

function NQCModels.potential!(model::AndersonHolsteinModel, V::Hermitian, r::AbstractMatrix)
    (;ħω, Ed, g) = model
    potential = ħω*first(r)^2/2
    V.data .= [potential sqrt(model.Γ/2π); sqrt(model.Γ/2π) potential + Ed + sqrt(2)*g*first(r)]
    return nothing
end

function NQCModels.derivative!(model::AndersonHolsteinModel, D::AbstractMatrix{<:Hermitian}, r::AbstractMatrix)
    (;ħω, g) = model
    D[1].data .= [ħω*first(r) 0; 0 ħω*first(r) + sqrt(2)*g]
    return nothing
end

model = AndersonHolsteinModel(ħω, Ed, g, Γ)

sim = Simulation{BCME}(atoms, model; temperature=kT, bandwidth=100.0)
```

The initial conditions are centred on the equilibrium position of the charged state ``x_1 = -\sqrt{2} g / \hbar\omega``:

```@example bcme
βω = ħω / kT
σ = sqrt(1 / βω)
r_dist = hcat(Normal(x1, σ))
v_dist = VelocityBoltzmann(kT, atoms.masses, (1,1))
distribution = DynamicalDistribution(v_dist, r_dist, (1,1)) * PureState(1, Diabatic())

output = run_dynamics(sim, (0.0, 5000.0), distribution;
    trajectories=50,
    output=OutputPosition,
    abstol=1e-10, reltol=1e-10,
    saveat=250
)
```

The mean position over time, normalised by ``x_1``, shows the equilibration of the system:

```@example bcme
using Statistics: mean
using Plots

times = output[1][:Time]
avg_pos = [mean(o[:OutputPosition][i][1] for o in output) / x1 for i in eachindex(output[1][:OutputPosition])]
plot(times, avg_pos; xlabel="t", ylabel="⟨x⟩ / x₁", label="BCME")
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
