```@setup logging
@info "Expanding src/NQCCalculators/overview.md..."
start_time = time()
```

# NQCCalculators.jl

`NQCCalculators` exists as the intermediary between `NQCModels` and `NQCDynamics`. 
It is home to all of the utility functions designed to calculate properties of single particle, 
multi-particle and ring-polymer systems interacting under the complete range of model potentials 
and electronic Hamiltonians implemented in `NQCModels`.
In addition, the module defines the data-structures, `caches`, that NQCDynamics interacts with for performing
simulation processes that depend on model quantities like the total energy derivative or the nonadiabtic 
coupling terms.
These `caches` are updated once per time-step during the dynamics simulation and are then subsequently used as a local
address in memory that contains all of the relavent quantities for a given dynamics method.

!!! note "What is the purpose of NQCCalculators?"
    `NQCCalculators.jl` formerly existed as a part of `NQCDynamics.jl` where it was concieved primarily as a way of
    calculating and caching relevant quantities for dynamics methods. While it continues to serve this purpose
    as a seperate module, it can now be more flexibly used in conjuction with `NQCModels.jl` as a collection of
    preimplemted methods for calculating the properties and observables of systems defined using the `NQCModels.jl` package.

# Using `Caches` - a simple example

In the overview of [`NQCModels`](@ref) we introduced a couple of simple [QuantumModels](@ref QuantumModels.QuantumModel)
and showed how to calculate some of their core properties with the `potential` and `derivative` functions. Here we will use the 
`model` object to build a `cache` that can calculate, manipulate and store more sophisitcated properties of the system.

```@example QuantumModel_Cache
using NQCCalculators
using NQCModels
using Unitful, UnitfulAtomic

#Creating a model
model = DoubleWell()

#Creating an empty cache with 3 atoms that propagate according to the double well potential
DoubleWell_Cache = Create_Cache(model, 3, Float64)
```

Here we have created an empty cache with all of the data-structures nessecary to contain quantities that define 3 independent atoms 
propagating under a double well potential. Notice how we didn't need to define the types of atoms in the system, only how many there are,
this is because the `cache` only stores data that relates to the electronic structure of the system and therefore depends only on a set
of atomic coordinates associated to positions within that electronic structure.

The `cache` itself is a `struct` that contains,
- the `model`
- the `potential` matrix
- the `potential` eigenvalues
- the `potential` eigenvectors
- the `derivative` matrices
- the `adiabatic derivative` matrices
- the `nonadiabtic coupling` matrices
- a temporary matrix for ad-hoc assignment during calculations

!!! info "Not all Caches store the same data"
    In general, the data the `cache` stores depends on the type of `model` and whether or not the atoms are being modelled with `ring polymers`.
    This is important for minimising memory overhead and customising the calculator routines to work most efficiently on each different
    system. More details about what each type of `cache` contains can be found in the relevant section of the `NQCCalculators` docs.

Now that we have created an empty `cache` to store all of the useful model data, we need to call calculation routines that will populate it with
the data we wish to store. This can be done easily by calling the `update_cache!` function which takes as input our `cache` and a set of nuclear 
coordintates and then populates the `cache` object we passed in with all of the values associated with that set of coordinates. 

For this example we will pass the function `DoubleWell_Cache` and a one-dimensional position value for each of the 3 atoms in our system.

```@example QuantumModel_Cache
#Populating DoubleWell_Cache with the values for a collection of atoms situated
#at 1.0, 2.0 and 3.0 angstroms
r = [austrip(1.0u"Å") austrip(2.0u"Å") austrip(3.0u"Å")]
update_cache!(DoubleWell_Cache, r)

#checking that update_cache!() has altered the values stored in DoubleWell_Cache
println(DoubleWell_Cache)

#Individual quantities can be extracted from the Cache by calling Cache.$quantity
ad_derivative = DoubleWell_Cache.adiabatic_derivative
eigvals = DoubleWell_Cache.eigen.values
eigvecs = DoubleWell_Cache.eigen.vectors
```

`update_cache!` is an easy and efficient way of calculating all of the relavent data associated to a model acting at a set of given coordinates. 
If, however, we are only interested in calculating a specific datapoint at a given positon, we can use the `update_$(quantity)!` functions.
Similarly to `update_cache!` they only require a `cache` and a position matrix, but these functions will only update the specific `$quantity`
for which they are named.

```@example QuantumModel_Cache
r_2 = [austrip(0.9u"Å") austrip(2.2u"Å") austrip(4.1u"Å")]
update_eigen!(DoubleWell_Cache, r_2)

#we can check that the new eigenvalues and eigenvectors are difference to before
eigvals = DoubleWell_Cache.eigen.values
eigvecs = DoubleWell_Cache.eigen.vectors

#we can also check that they are the only values that have been updated
ad_derivative = DoubleWell_Cache.adiabatic_derivative
```

Alternatively, the `evaluate_$(quantity)` functions can be used to calculate a specific datapoint at a given positon, without updating the cache. The "evaluated" `$quantity` for which the function is named is output to an assigned variable, leaving the cache unchanged. Again, these functions only require a `cache` and a position matrix.

```@example QuantumModel_Cache
r_3 = [austrip(0.8u"Å") austrip(2.4u"Å") austrip(4.2u"Å")]
DoubleWell_eigen = evaluate_eigen(DoubleWell_Cache, r_3)

#we can check that the new evaluated eigenvalues and eigenvectors are different to those currently present in the cache
eigvals_eval = DoubleWell_eigen.values
eigvals_cache = DoubleWell_Cache.eigen.values

eigvecs_cache = DoubleWell_eigen.vectors
eigvecs_cache = DoubleWell_Cache.eigen.vectors
```

!!! info "update!" vs evaluate" functions
    Notice that the `update_$(quantity)!` functions have an exclimation point, `!`, whereas the `evaluate_$(quantity)` functions do not. 
    This is common Julia notation for indicating an **in-place function**. An in-place function does not return a variable but instead updates one of the provided arguments "in-place". 
    As such, the `update_$(quantity)!` functions are written as inplace functions as they "update" fields in the cache that they are provided with, without returning a new variable. For the `evaluate_$(quantity)` functions however, we want a new variable to be returned and for the cache in its original state to be preserved, so regular (non-in-place) functions are used (indicated by no `!` in the name) that return new variables containing the "evaluated" quantity. 

