```@setup logging
@info "Expanding src/api/NQCModels/adiabaticmodels.md..."
start_time = time()
```

# AdiabaticModels

```@autodocs
Modules=[NQCModels.AdiabaticModels]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
