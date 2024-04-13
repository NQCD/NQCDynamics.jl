```@setup logging
@info "Expanding src/api/NQCDynamics/structure.md..."
start_time = time()
```
# Structure
This submodule contains utility functions to analyse and modify atomic structure data, such as interatomic distances, centres of mass, both with and without support for periodic boundary conditions. 

These functions can be used to build more sophisticated output functions, or for basic analysis of simulation results in post. 

This module **doesn't contain**:
- Basic definitions of atomic structures (e.g. `Atoms`, `PeriodicCell`, ...). These are defined in [`NQCBase`](@ref).
- Functions to generate atomic structures. These should be added to [`NQCDynamics.InitialConditions`](@ref).
- Analysis functions for specific systems (e.g. molecules on surfaces). These should be added to [`NQCDynamics.Analysis`](@ref). 


## Method reference


```@autodocs
Modules=[NQCDynamics.Structure]
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
