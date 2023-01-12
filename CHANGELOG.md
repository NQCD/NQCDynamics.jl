# NQCDynamics.jl changelog

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
