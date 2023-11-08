```@setup logging
@info "Expanding src/api/NQCDynamics/estimators.md..."
start_time = time()
```

# Estimators

```@autodocs
Modules=[NQCDynamics.Estimators]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
