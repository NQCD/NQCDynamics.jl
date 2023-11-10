```@setup logging
@info "Expanding src/api/NQCDynamics/dynamicsutils.md..."
start_time = time()
```

# DynamicsUtils

```@autodocs
Modules=[NQCDynamics.DynamicsUtils]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
