```@setup logging
@info "Expanding src/api/NQCDynamics/dynamicsmethods.md..."
start_time = time()
```

# DynamicsMethods

```@autodocs
Modules=[NQCDynamics.DynamicsMethods]
```

## ClassicalMethods

```@autodocs
Modules=[NQCDynamics.DynamicsMethods.ClassicalMethods]
```

## MappingVariableMethods

```@autodocs
Modules=[NQCDynamics.DynamicsMethods.MappingVariableMethods]
```

## SurfaceHoppingMethods

```@autodocs
Modules=[NQCDynamics.DynamicsMethods.SurfaceHoppingMethods]
```

## EhrenfestMethods

```@autodocs
Modules=[NQCDynamics.DynamicsMethods.EhrenfestMethods]
```

## IntegrationAlgorithms

```@autodocs
Modules=[NQCDynamics.DynamicsMethods.IntegrationAlgorithms]
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
