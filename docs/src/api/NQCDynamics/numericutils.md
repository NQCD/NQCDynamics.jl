```@setup logging
@info "Expanding src/api/NQCDynamics/numericutils.md..."
start_time = time()
```

# Numerical utilities

```@autodocs
Modules=[NQCDynamics.FastDeterminant]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
