```@setup logging
@info "Expanding src/api/NQCModels/adiabaticmodels.md..."
start_time = time()
```

# ClassicalModels

```@autodocs
Modules=[NQCModels.ClassicalModels]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
