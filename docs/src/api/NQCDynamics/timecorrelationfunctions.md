```@setup logging
@info "Expanding src/api/NQCDynamics/timecorrelationfunctions.md..."
start_time = time()
```

# TimeCorrelationFunctions

```@autodocs
Modules=[NQCDynamics.TimeCorrelationFunctions]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
