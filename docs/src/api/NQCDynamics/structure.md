```@setup logging
@info "Expanding src/api/NQCDynamics/structure.md..."
start_time = time()
```
# Structure

```@autodocs
Modules=[NQCDynamics.Structure]
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
