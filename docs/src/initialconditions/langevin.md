```@setup logging
@info "Expanding src/initialconditions/langevin.md..."
start_time = time()
```
# [Thermal Langevin dynamics](@id langevin-sampling)

Thermal initial conditions can be obtained directly from a Langevin dynamics simulations.
See the [Langevin dynamics page](@ref langevin-dynamics) for more info.
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
