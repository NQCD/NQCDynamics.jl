```@setup logging
@info "Expanding src/api/NQCBase/nqcbase.md..."
start_time = time()
```

# NQCBase

```@autodocs
Modules=[NQCBase]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
