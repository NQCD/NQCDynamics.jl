```@setup logging
@info "Expanding src/api/NQCDynamics/analysis.md..."
start_time = time()
```
# Analysis

```@autodocs
Modules=[NQCDynamics.Analysis]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```