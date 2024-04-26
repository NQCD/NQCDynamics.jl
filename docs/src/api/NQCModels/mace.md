```@setup logging
@info "Expanding src/api/NQCModels/mace.md..."
start_time = time()
```

# MACE interface

## Possible output data

```@autodocs
Modules=[NQCModels.AdiabaticModels.MACE]
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
