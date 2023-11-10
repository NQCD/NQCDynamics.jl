```@setup logging
@info "Expanding src/api/NQCDynamics/ringpolymers.md..."
start_time = time()
```

# RingPolymers

```@autodocs
Modules=[NQCDynamics.RingPolymers]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
