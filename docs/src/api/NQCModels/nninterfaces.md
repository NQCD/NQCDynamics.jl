```@setup logging
@info "Expanding src/api/NQCModels/nninterfaces.md..."
start_time = time()
```

# NNInterfaces

```@autodocs
Modules=[NNInterfaces]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
