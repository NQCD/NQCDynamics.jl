```@setup logging
@info "Expanding src/api/NQCDInterfASE/ase-interface.md..."
start_time = time()
```

# NQCDInterfASE

```@autodocs
Modules=[NQCDInterfASE]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
