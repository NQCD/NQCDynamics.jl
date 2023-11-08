```@setup logging
@info "Expanding src/api/NQCDynamics/nonadiabaticmoleculardynamics.md..."
start_time = time()
```
# NQCDynamics

```@autodocs
Modules=[NQCDynamics]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
