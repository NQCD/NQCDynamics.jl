# Implementing a new model

[NonadiabaticModels.jl](@ref) aims to provide a unified interface for defining both model
Hamiltonians for nonadiabatic dynamics simulations.

Here we shall walk through the implementation of a few differnet types of model to
hopefully show the simplicity of the interface.

## Abstract types

Julia's [abstract type](https://docs.julialang.org/en/v1/manual/types/#man-abstract-types)
system can be likened to the inheritance concept from object-oriented programming or the
trait system from [Rust](https://doc.rust-lang.org/book/ch10-02-traits.html).
It allows us to defined shared behaviour for a groups of structs and allows us to define
a common set of functions that all of the concrete types must implement.

Here, the top level abstract type is the [`Model`](@ref NonadiabaticModels.Model),
under which all of our models must fall.
The second tier below this includes the two abstract types
[`AdiabaticModel`](@ref NonadiabaticModels.AdiabaticModels) and
[`DiabaticModel`](@ref NonadiabaticModels.DiabaticModels).
These form the two distinct branches within the NonadiabaticModels type hierachy and the
shared behaviour across the branches is minimal.
The [`AdiabaticModel`](@ref NonadiabaticModels.AdiabaticModels.AdiabaticModel)
describes the familiar
form from molecular dynamics that provides
a single provides a single potential energy surface.
The [`DiabaticModel`](@ref NonadiabaticModels.DiabaticModels.DiabaticModel) instead provides 
multiple potential energy surfaces with couplings between them, as implied by the name,
these are in the diabatic representation.
If the desired model does not fall under either of these branches, a new abstract type
should be created.

### Minor branches

Under the two main branches there also exists some specialised abstract types that are
useful in some cases, such as when using many electronic states, or when implementing
extra functions is required. See the docstrings for more info:

- [`AdiabaticFrictionModel`](@ref NonadiabaticModels.FrictionModels.AdiabaticFrictionModel)
- [`LargeDiabaticModel`](@ref NonadiabaticModels.DiabaticModels.LargeDiabaticModel)
- [`DiabaticFrictionModel`](@ref NonadiabaticModels.DiabaticModels.DiabaticFrictionModel)

## Example implementations

With the information above it should be possible to select the appropriate abstract type
for your new model.
To implement a new model for your chosen abstract type, you should first take a look
at the docstring for the abstract type (click on the links above).
This will include a list of the functions that should be implemented along with an example 
implementation.

If this is not satisfactory, take a look into the source code of the `NonadiabaticModels.jl`
package to see how the analytic models have been implemented. 
If you have any issues or questions about implementing a new model, open up an issue on
Github and we can work together to resolve the problem. 
