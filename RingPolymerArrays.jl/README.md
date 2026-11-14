# RingPolymerArrays

[![CI](https://github.com/NQCD/RingPolymerArrays.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/NQCD/RingPolymerArrays.jl/actions/workflows/CI.yml)

RingPolymerArrays.jl provides the `RingPolymerArray{T}`, an `AbstractArray{T,3}`, useful for representing ring polymer systems using the imaginary-time path integral formalism.

A `RingPolymerArray` can be used to store positions or momenta for ring polymer systems, where the three dimensions represent each degree of freedom per atom, each atom and each ring polymer bead.

For example, this `RingPolymerArray` represents a 3D system with 2 atoms, and 4 ring polymer beads:
```julia
julia> RingPolymerArray(reshape(1:24, 3, 2, 4))
3×2×4 RingPolymerArray{Int64}:
[:, :, 1] =
 1  4
 2  5
 3  6

[:, :, 2] =
 7  10
 8  11
 9  12

[:, :, 3] =
 13  16
 14  17
 15  18

[:, :, 4] =
 19  22
 20  23
 21  24
```

In a system with disparate particle masses, it is desirable to treat heavier atoms classically, without any replicas.
This is made possible by the `classical` keyword argument in the constructor that allows us to specify which atoms should be treated classically.
Here, the first atom will have a fixed value for all beads:
```julia
julia> a = RingPolymerArray{Float64}(reshape(1:24, 3, 2, 4); classical=[1])
3×2×4 RingPolymerArray{Float64}:
[:, :, 1] =
 1.0  4.0
 2.0  5.0
 3.0  6.0

[:, :, 2] =
 1.0  10.0
 2.0  11.0
 3.0  12.0

[:, :, 3] =
 1.0  16.0
 2.0  17.0
 3.0  18.0

[:, :, 4] =
 1.0  22.0
 2.0  23.0
 3.0  24.0
```
These values will remain fixed, even after attempted modification:
```julia
julia> a .+= rand(3, 2, 4)
3×2×4 RingPolymerArray{Float64}:
[:, :, 1] =
 1.78368  4.85881
 2.90528  5.81841
 3.57008  6.31728

[:, :, 2] =
 1.78368  10.6239
 2.90528  11.9604
 3.57008  12.4336

[:, :, 3] =
 1.78368  16.1968
 2.90528  17.8464
 3.57008  18.8612

[:, :, 4] =
 1.78368  22.7994
 2.90528  23.1735
 3.57008  24.7504
 ```
 This works by internally storing only a single copy of the classical positions and allowing `setindex!` to operate only on the first bead for atoms labelled as classical.
 
## Normal mode transformation

Alongside the `RingPolymerArray` we also provide the `NormalModeTransformation{T}` which can be used to convert any `AbstractArray{T,3}` of ring polymer coordinates to and from normal mode coordinates:
```julia
julia> transform = NormalModeTransformation{Float64}(4);

julia> transform_to_normal_modes!(fill(1.0, (2,3,4)), transform)
2×3×4 Array{Float64, 3}:
[:, :, 1] =
 2.0  2.0  2.0
 2.0  2.0  2.0

[:, :, 2] =
 -8.65956e-17  -8.65956e-17  -8.65956e-17
 -8.65956e-17  -8.65956e-17  -8.65956e-17

[:, :, 3] =
 0.0  0.0  0.0
 0.0  0.0  0.0

[:, :, 4] =
 2.59787e-16  2.59787e-16  2.59787e-16
 2.59787e-16  2.59787e-16  2.59787e-16
 
julia> transform_from_normal_modes!(fill(1.0, (2,3,4)), transform)
2×3×4 Array{Float64, 3}:
[:, :, 1] =
 1.70711  1.70711  1.70711
 1.70711  1.70711  1.70711

[:, :, 2] =
 -0.707107  -0.707107  -0.707107
 -0.707107  -0.707107  -0.707107

[:, :, 3] =
 0.292893  0.292893  0.292893
 0.292893  0.292893  0.292893

[:, :, 4] =
 0.707107  0.707107  0.707107
 0.707107  0.707107  0.707107
```
When using a `RingPolymerArray`, the normal mode transformation is applied to both classical
and quantum atoms, where the classical atoms will have the value of the scaled centroid.
In normal mode coordinates, the centroid is given by the first bead coordinate divided by `sqrt(nbeads)`.
