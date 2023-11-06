```@setup logging
@info "Expanding src/api/NQCDistributions/nqcdistributions.md..."
start_time = time()
```

# NQCDistributions

```@autodocs
Modules=[NQCDistributions]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
