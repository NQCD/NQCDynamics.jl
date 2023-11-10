```@setup logging
@info "Expanding src/api/NQCDynamics/dynamicsoutputs.md..."
start_time = time()
```

# DynamicsOutputs

Here are all the functions that you can specify in the output tuple when using
`run_dynamics`.
To add more, simply add a new function in the `DynamicsOutputs` module. 
```@autodocs
Modules=[NQCDynamics.DynamicsOutputs]
Private=false
```

## Internals

```@autodocs
Modules=[NQCDynamics.DynamicsOutputs]
Public=false
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
