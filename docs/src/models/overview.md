# Overview

```@docs
Models
```

## Abstract types
```@docs
Models.Model
Models.AdiabaticModel
Models.DiabaticModel
Models.FrictionModel
```

## Functions to implement
To create a new model, depending on the type, it should implement some of these functions.
```@docs
Models.potential!
Models.derivative!
Models.friction!
```