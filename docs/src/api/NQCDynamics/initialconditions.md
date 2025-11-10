```@setup logging
@info "Expanding src/api/NQCDynamics/initialconditions.md..."
start_time = time()
```

# InitialConditions

```@autodocs
Modules=[NQCDynamics.InitialConditions]
```

## ThermalMonteCarlo

```@autodocs
Modules=[NQCDynamics.InitialConditions.ThermalMonteCarlo]
```

## QuantisedDiatomic

```@autodocs
Modules=[NQCDynamics.InitialConditions.QuantisedDiatomic]
```

## MetropolisHastings

```@autodocs
Modules=[NQCDynamics.InitialConditions.MetropolisHastings]
```

## ConfigureAtomic

```@autodocs
Modules=[NQCDynamics.InitialConditions.ConfigureAtomic]
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
