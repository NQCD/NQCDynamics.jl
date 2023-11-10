```@setup logging
@info "Expanding src/api/NQCDynamics/calculators.md..."
start_time = time()
```
# Calculators

```@autodocs
Modules=[NQCDynamics.Calculators]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
