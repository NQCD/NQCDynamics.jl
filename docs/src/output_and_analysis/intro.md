# [Simulations outputs and analysis functions](@id output_and_analysis)

NQCDynamics.jl is able to output a variety of common observables when running simulations with the `run_dynamics` function. Further output functions can be implemented easily as well. 

In this section you will find an overview of all available output types, as well as explanations of some common analysis functions in the realm of surface chemistry which we have implemented in the package. 

## DynamicsOutputs

In many examples within this documentation, you have seen calls to `run_dynamics`:

```julia
ensemble = run_dynamics(sim, tspan, distribution;selection=1:20,
    dt=0.1u"fs", output=OutputPosition, trajectories=20, callback=terminate)
```

Within `run_dynamics`, the `output` argument specifies the desired output values for each trajectory. `output` can either be given as a single function, or as a tuple of multiple output functions, for example:

```julia
output=OutputPosition # or
output=(OutputPosition, OutputVelocity, OutputKineticEnergy)

ensemble[3][:OutputPosition] # will output the positions at all timesteps in trajectory 3
```

Every output type is a function which can use the [`DynamicsVariables`](@ref) and [`Simulation`](@ref) values of the respective trajectory, allowing you to create custom output types of any kind. See the [developer documentation] for more information on how to implement a custom output type. 

You can find an overview of all available output types in the [`DynamicsOutputs`](@ref NQCDynamics.DynamicsOutputs) API. 

## Analysis functions

The Analysis submodule in NQCDynamics contains functions commonly used in the analysis of trajectories to make the analysis of existing trajectories easier. 
Ideally, most observable quantities could be implemented with a combination of [`DynamicsOutputs`](@ref NQCDynamics.DynamicsOutputs) and `Reduction` types, however we might want to data from existing ensemble simulations where re-running the entire set of trajectories is impractical. 

As a result, most functions in the `Analysis` submodule are also implemented as a `DynamicsOutput`. 

### Convenient functions for periodic structures
[NQCDynamics.Structure](@ref NQCDynamics.Structure) contains several useful functions for periodic structures, such as `pbc_distance, pbc_center_of_mass`. 

These functions take into account periodic copies of the atoms in question, returning the respective values for the closest set of periodic copies. 

### Analysis of diatomic molecules
[NQCDynamics.Analysis.Diatomic](@ref NQCDynamics.Analysis.Diatomic) 
