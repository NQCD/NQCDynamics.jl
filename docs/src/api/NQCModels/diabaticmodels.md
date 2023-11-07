```@setup logging
@info "Expanding src/api/NQCModels/diabaticmodels.md..."
start_time = time()
```

# DiabaticModels

```@autodocs
Modules=[NQCModels.DiabaticModels]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
