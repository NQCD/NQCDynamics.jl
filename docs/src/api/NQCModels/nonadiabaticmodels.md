```@setup logging
@info "Expanding src/api/NQCModels/nonadiabaticmodels.md..."
start_time = time()
```

# NQCModels

```@autodocs
Modules=[NQCModels]
```



```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
