```@setup logging
@info "Expanding src/api/NQCDynamics/calculators.md..."
start_time = time()
```
# NQCCalculators

```@autodocs
Modules=[NQCCalculators]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
