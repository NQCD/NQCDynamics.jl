```@setup logging
@info "Expanding src/NQCModels/combining_models.md..."
start_time = time()
```

# Composing multiple models

## Why?

For some non-adiabatic dynamics methods, we would like to calculate certain quantities using a number of different approximations or different machine learning model interfaces. 
However, a model for a potential energy surface does not and should not be able to provide e.g. an electronic friction tensor. 
Instead, we can use the modular nature of NQCModels.jl to combine different components to a non-adiabatic dynamics method into a single model that can be used to set up a [`Simulation`](@ref)

For example, in order to run [Molecular Dynamics with Electronic friction](@ref mdef-dynamics), we need to use a [`Model`](@ref) which supports the following methods:
- `potential` to determine the total energy
- `derivative` to determine interatomic forces
- `friction` to determine the electronic friction tensor

However, most models for potential energy surfaces only provide `potential` and `derivative`, whereas `friction!` is provided by an `ElectronicFrictionProvider`. 

To combine different models together for a system, we can use [`Subsystem`](@ref)s to define which atoms in the system a given model should be applied to, and a [`CompositeModel`](@ref) to combine the different models into one. 
If different [`Subsystem`](@ref)s should experience different effective temperatures, a [`Thermostat`](@ref) can be provided when initialising a [`Simulation`](@ref). 

![An overview of the different components used to combine models in NQCModels.jl](../assets/compositemodels/struct-explainer.svg)

### Example: Combining an `AdiabaticModel` with `ElectronicFrictionProvider`s

```@example compositemodels
using NQCDynamics

atoms = vcat([:H, :H], [:Cu for _ in 1:54]) # create example surface structure
positions = rand(3,56) # Placeholder for atomic positions
pes_model = Subsystem(Free(3)) # PES model with potential and derivative
friction_model = Subsystem(RandomFriction(3)) # Friction model supporting friction! only

println("PES model: \n", pes_model, "\nFriction model:", "friction_model")

complete_model = CompositeModel(pes_model, friction_model)
```

As shown above, we have combined the `RandomFriction` provider with a model potential to give a total model which can be used in a simulation. 

```@example compositemodels
println("Potential: ", NQCModels.potential(complete_model, positions))
println("Derivative: ", NQCModels.derivative(complete_model, positions))
println("Friction: ", NQCModels.friction(complete_model, positions))
```
