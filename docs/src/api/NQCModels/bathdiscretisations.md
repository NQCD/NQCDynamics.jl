```@setup logging
@info "Expanding src/api/NQCModels/adiabaticmodels.md..."
start_time = time()
```

# BathDiscretisations

```@autodocs
Modules=[NQCModels.BathDiscretisations]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
