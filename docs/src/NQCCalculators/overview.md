```@setup logging
@info "Expanding src/NQCCalculators/overview.md..."
start_time = time()
```

# NQCCalculators.jl

`NQCCalculators` exists as the intermediary between `NQCModels` and `NQCDynamics`. 
It is home to all of the utility functions designed to calculate properties of single particle, 
multi-particle and ring-polymer systems interacting under the complete range of model potentials 
and electronic Hamiltonians implemented in `NQCModels`.
In addition, the module defines the data-structures, `Caches`, that NQCDynamics interacts with for performing
simulation processes that depend on model quantities like the total energy derivative or the non-adiabtic 
coupling terms.
These `Caches` are updated once per time-step in the dynamics simulation and are then subsequently used as a local
address in memory that contains all of the relavent quantities for a given dynamics method.

!!! note

  `NQCCalculators.jl` formerly existed as a part of `NQCDynamics.jl` where it was concieved primarily as a way of 
  calculating and caching relevant quantities for dynamics methods. While it continues to serve this purpose 
  as a seperate module, it can now be more flexibly used in conjuction with `NQCModels.jl` as a collection of 
  preimplemted methods for calculating the properties and observables of instances of the systems defined
  within the `NQCModels.jl` package.

# Using `Caches` - a simple example

In the [overview of `NQCModels`](@ref) we introduced a couple of simple [QuantumModels](@ref NQCModels.QuantumModels.QuantumModel)
and showed how to calculate some of their core properties with the `potential` and `derivative` functions. Here we will use the 
`model` object to build a `Cache` which will be used to calculate, manipulate and store more sophisitcated properties of the system.

```@example QuantumModel_Cache
using NQCCalculators
using NQCModels
using Unitful, UnitfulAtomic

#Creating a model
model = DoubleWell()

#Creating an empty Cache with the new model and 3 atoms to evolve under it
DoubleWell_Cache = Create_Cache(model, 3, Float64)

#Populating the new Cache with the values for atoms 1.0, 2.0 and 3.0 angstroms from the surface
update_cache!(DoubleWell_Cache, [austrip(1.0u"Å") austrip(2.0u"Å") austrip(3.0u"Å")])
```
