```@setup logging
@info "Expanding src/api/NQCDynamics/ensembles.md..."
start_time = time()
```

# Ensembles

```@autodocs
Modules=[NQCDynamics.Ensembles]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
