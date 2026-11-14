# NQCModels.jl

[![CI](https://github.com/NQCD/NQCModels.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/NQCD/NQCModels.jl/actions/workflows/CI.yml)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nqcd.github.io/NQCDynamics.jl/dev/NQCModels/overview/)

This package provides an interface for defining models for nonadiabatic dynamics.
Primarily, the package is intended to be used with [NQCDynamics.jl](https://github.com/NQCD/NQCDynamics.jl/)
but it is equally possible to use any of the models separately.

For details on the included models and implementation take a look at the [documentation](https://nqcd.github.io/NQCDynamics.jl/dev/NQCModels/overview/).

## Features

- **Unified interface** for classical and quantum models
- **Extensive model library** including Tully models, Ananth models, and many others
- **GPU acceleration** via CUDA.jl (NVIDIA) and Metal.jl (Apple Silicon) for many models (requires Julia 1.9+)
- **Composable models** for system-bath dynamics and model combinations
- **Full-dimensional interfaces** including machine learning models and AtomsCalculators

See the [GPU Acceleration](https://nqcd.github.io/NQCDynamics.jl/dev/NQCModels/gpu_support/) documentation for details on GPU support.

> [!NOTE]
> [ACEregistry](https://github.com/ACEsuit/ACEregistry) must be installed if you are running Julia < v1.11. 
