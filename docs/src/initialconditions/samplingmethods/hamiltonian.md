# Thermal Hamiltonian Monte Carlo

Our implementation of Hamiltonian Monte Carlo (HMC) is a light wrapper around the
[AdvancedHMC.jl](https://github.com/TuringLang/AdvancedHMC.jl) package.
If you want to learn about the HMC theory, refer to the references that they provide.

Currently this works only for systems with classical nuclei
(i.e. `Simulation` not `RingPolymerSimulation`)
but it should be possible to extend it in the future.

## Example

```@example hmc
using NonadiabaticMolecularDynamics
using Unitful
using UnitfulAtomic

sim = Simulation(Atoms([:H, :H, :C, :C]), Harmonic(dofs=3); temperature=300u"K")
r0 = randn(size(sim))
chain, stats = InitialConditions.ThermalMonteCarlo.run_advancedhmc_sampling(sim, r0, 1e4)
nothing # hide
```

Here we should see that the energy expectation for the generated ensemble matches
with the equipartion theorem:
```@repl hmc
Estimators.@estimate potential_energy(sim, chain)
austrip(sim.temperature) * 3 * 4 / 2
```
