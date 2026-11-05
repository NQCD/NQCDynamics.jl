# NQCDynamics.jl changelog

## Yet to be named GPU version of NQCD

### NQCDynamics.jl
- Removed `TemperatureSetting` and replaced with the ability to have a single temperature function return a length(Atoms) Vector for the temperature. 

### NQCModels.jl
- Added experimental Metal.jl and CUDA.jl support for a number of models. 
- ![bugfix][badge-bugfix] `ClassicalModel`s are now consistently multi-particle capable: `potential(model, R)` always returns a scalar summed over every atom (column of `R`), and `derivative`/`derivative!` always return/fill a `dimensions * atoms` array of scalars, with a generic default provided for models that treat every atom independently.
- ![bugfix][badge-bugfix] `QuantumModel`s are now consistently multi-particle capable: `potential(model, R)` always returns a single `nstates * nstates` Hermitian matrix summed over every independent atom, and `derivative`/`derivative!` always return/fill a `dimensions * atoms` array of Hermitian matrices. A generic default `potential!` (mirroring the existing `derivative!` default) was added so 1-dof-per-atom models only need to implement a single-atom `potential!(model, V, r::Number)`.

  The following models were found to violate this contract and have been fixed:

  | Model | Issue |
  |---|---|
  | `TullyModelOne`/`Two`/`Three` | `potential(model, R)` always returned an all-zero matrix; no matrix-level `potential!` existed to be reached by the (now added) generic default. |
  | `Scattering1D` | Same missing matrix-level entry point as Tully, plus `potential!` double-counted the `V₀` contribution at the molecule-state diagonal entry, inconsistent with its own non-mutating `potential`. |
  | `ThreeStateMorse` | `potential` never summed over atoms (only read the first atom's position); its own non-mutating `derivative`/`derivative!` returned/filled a bare `Hermitian` instead of the required `Matrix{Hermitian}` of shape `ndofs * natoms`. |
  | `MiaoSubotnik` | Same atom-summing and derivative-shape bugs as `ThreeStateMorse`. |
  | `OuyangModelOne` | Same atom-summing and derivative-shape bugs as `ThreeStateMorse`, plus `potential!` was missing the diagonal-states-fill line present in the non-mutating `potential`, so the in-place version was incomplete even for a single atom. |
  | `GatesHollowayElbow` | `potential`/`potential!`/`derivative!` only ever read the first atom's coordinates; extra atoms were silently ignored instead of summed. |
  | `ErpenbeckThoss` | `potential`'s off-diagonal/second-state terms were not summed over atoms (only the Morse-derived term was); `derivative!` had an out-of-bounds indexing bug (`D[i,i]` instead of `D[1,i]`) and passed a `Vector` where an `AbstractMatrix` was required, causing a `MethodError` even for a single atom. |
  | `AndersonHaldane` | Called `ClassicalModels.Morse`/`Logistic` with a raw scalar, which their API never accepted, causing a `MethodError` regardless of atom count. |
  | `DoubleWell` | `derivative!` indexed per-atom parameter vectors (`mass`, `ω`, `γ`) directly instead of broadcasting, causing a `BoundsError` for more atoms than parameters provided; separately, the electronic coupling term `Δ` did not scale with atom count in `potential`/`potential!`, unlike every other term. |
  | `AndersonHolstein` | Cached an `impurity_derivative` buffer sized for exactly 1 atom, causing a crash for any other atom count; also reused a persistent `impurity_potential` buffer across calls without zeroing it, causing incorrect accumulation once the wrapped model's `potential!` correctly accumulates (`+=`) into `V`. |
  | `BosonBath` | `potential!` reassigned the local `V` argument instead of mutating it in place (a silent no-op). |

  `AnanthModelOne`/`AnanthModelTwo`, `SpinBoson`, `WideBandBath` and `AdiabaticStateSelector` were already correctly multi-particle capable and needed no changes.
- ![bugfix][badge-bugfix] Fixed a method ambiguity in `CompositeFrictionModel.derivative!` introduced by the `ClassicalModel` generic default above.

## Version `v1.0.3`
- ![Maintenance][badge-maintenance] Updated documentation to include JLD2.jl installation instructions [#457][github-457]
- Fixed missing integrator allocation for FSSH, which lead to incorrect simulation results. 

### NQCModels.jl v1.1.0
- `LuHertlMaurerModel` renamed to `AndersonHaldane`

### NQCBase v1.1.0
Added `Structure` type which contains atoms, positions, cell and additional information, similar to other atomic structure manipulation frameworks. 

This lays the groundwork to multiple dispatch off of atomic structure information, which could be used to provide plotting recipes to Makie.jl / Plots.jl in future. [Related PR](https://github.com/NQCD/NQCBase.jl/pull/30)

### NQCDInterfASE v1.1.0
> [!IMPORTANT]
>  `convert_from_ase_atoms` now outputs a `Structure`. Previously, it would output a `Tuple` of `(atoms, positions, cell)`

### MACEModels v1.5.0
- Added support for newest version of MACE.
- Allow loading foundation models, multi-head models and mace-freeze models.
- Fixed issue where having multiple models with different dtypes loaded would cause errors depending on which model was created first.

### NQCCalculators v1.0.10
- Switched to FastLapackInterface

## Version `v0.13.0`

- ![enhancement][badge-enhancement] Added BCME method [#290][github-290]
- ![enhancement][badge-enhancement] Added CME method [#290][github-290]
- ![enhancement][badge-enhancement] Added Ehrenfest for Newns-Anderson models [#290][github-290]
- ![BREAKING][badge-breaking] New default algorithm for Ehrenfest [#290][github-290]
- ![Fix][badge-bugfix] Improved probability estimate for IESH [#290][github-290]
- ![Maintenance][badge-maintenance] Removed DiabaticIESH [#290][github-290]

## Version `v0.12.1`

- ![enhancement][badge-enhancement] Added Independent Electron Surface Hopping Method [#287][github-287]

## Version `v0.12.0`

- ![BREAKING][badge-breaking] Merged `run_trajectory` and `run_ensemble` into `run_dynamics` [#286][github-286]
- ![BREAKING][badge-enhancement] Output from `run_dynamics` is a `Dictionary` from `Dictionaries.jl` [#286][github-286]
- ![enhancement][badge-enhancement] `FileReduction(filename)` allows for output to be stored in HDF5 format [#286][github-286]

## Version `v0.11.0`

- ![Maintenance][badge-maintenance] Moved `NonadiabaticDistributions` to [NQCDistributions.jl]() [#275][github-275]
- ![BREAKING][badge-breaking] `BoltzmannVelocityDistribution` renamed to VelocityBoltzmann [#275][github-275]
- ![BREAKING][badge-breaking] `SingleState`, `ElectronicPopulation` renamed to `PureState`, `MixedState` [#275][github-275]
- ![enhancement][badge-enhancement] Created changelog!

[github-457]: https://github.com/NQCD/NQCDynamics.jl/pull/457
[github-290]: https://github.com/NQCD/NQCDynamics.jl/pull/290
[github-287]: https://github.com/NQCD/NQCDynamics.jl/pull/287
[github-286]: https://github.com/NQCD/NQCDynamics.jl/pull/286
[github-275]: https://github.com/NQCD/NQCDynamics.jl/pull/275

[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/feature-green.svg
[badge-enhancement]: https://img.shields.io/badge/enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/bugfix-purple.svg
[badge-security]: https://img.shields.io/badge/security-black.svg
[badge-experimental]: https://img.shields.io/badge/experimental-lightgrey.svg
[badge-maintenance]: https://img.shields.io/badge/maintenance-gray.svg
