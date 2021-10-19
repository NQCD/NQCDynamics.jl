# [Electronic friction models](@id models-friction)

To perform [molecular dynamics with electronic friction (MDEF)](@ref mdef-dynamics)
a specific type of model must be used,
these provide the friction tensor used to propagate the dynamics.

Since the electronic friction concept relies directly upon an approximate interaction
of the nuclei with the metal electrons, it is usually necessary to interface with *ab initio*
data.
As detailed in the [MDEF page](@ref mdef-dynamics), there are two ways to obtain friction
values, either from the local density friction approximation (LDFA), or from time-dependent
perturbation theory (TDPT).

## Analytic models

Since *ab initio* friction calculations are often expensive it is useful to investigate
have some models which we can use to test friction methods.
The [`DiabaticFrictionModel`](@ref NonadiabaticModels.DiabaticModels.DiabaticFrictionModel)
is provided such that the electronic friction can be directly evaluated for diabatic models.
These calculate the friction directly from the nonadiabatic couplings for the model.

## [CubeLDFAModel.jl](@id models-cubeldfa)

Our LDFA implementation is given in
[CubeLDFAModel.jl](https://github.com/NQCD/CubeLDFAModel.jl)
which takes a `.cube` file containing the electron density and evaluates the friction based
upon this local density.

The model works by fitting the LDA data provided by [Gerrits2020](@cite) which provides
the LDFA friction coefficient as a function of the Wigner-Seitz radius.
Through the following relation:
```math
r_s(\rho) = (\frac{3}{4\pi \rho (\mathbf{r_{i}})})^{1/3}
```
we can immediately connect the electron density to the radius, and then to the friction values.
Now it is a simple matter of evaluating the electron density at a given point and converting
the density into the friction value.
Optimally, this would be done via an *ab initio* calculation to get the electron density,
but this model instead uses a pre-computed `.cube` file to get the density with minimal cost.
This makes the assumption that the density does not change throughout the dynamics, or that
the surface is assumed to be frozen in place.

This graph shows how we interpolate the LDA data and evaluate the friction coefficient
as a function of the Wigner-Seitz radius.
![ldfa graph](../assets/figures/ldfa_graph.png)

The [reactive scattering example](@ref example-h2scattering) uses this model to investigate
the scattering of a diatomic molecule from a metal surface.

## NNInterfaces.jl

Another way to perform MDEF simulations is the use one of the models from
[`NNInterfaces.jl`](https://github.com/NQCD/NNInterfaces.jl/) which uses a neural network
to obtain the time-dependent perturbation theory friction from the atomic positions.
As with LDFA, one of these models is used in the
[reactive scattering example](@ref example-h2scattering).

