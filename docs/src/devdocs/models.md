```@setup logging
@info "Expanding src/devdocs/models.md..."
start_time = time()
```
# [Implementing a new model](@id devdocs-model)

[NQCModels.jl](@ref) aims to provide a unified interface for defining model
Hamiltonians for nonadiabatic dynamics simulations.

Here, we walk through the implementation of a few different types of model to
explain the interface.

## Abstract types

Julia's [abstract type](https://docs.julialang.org/en/v1/manual/types/#man-abstract-types)
system can be likened to the inheritance concept from object-oriented programming or the
trait system from [Rust](https://doc.rust-lang.org/book/ch10-02-traits.html).
It allows us to defined shared behaviour for a groups of `struct`s and allows us to define
a common set of functions that all of the concrete types must implement.

In [NQCModels.jl](@ref) the top level abstract type is the [`Model`](@ref NQCModels.Model),
under which all of our models must fall.
The second tier below this includes the two abstract types
[`ClassicalModel`](@ref ClassicalModels) and
[`QuantumModel`](@ref QuantumModels).
These form the two distinct branches within the NQCModels type hierarchy and the
shared behaviour across the branches is minimal.
The [`ClassicalModel`](@ref ClassicalModels.ClassicalModel)
describes the familiar form from molecular dynamics that provides a single potential energy surface.
The [`QuantumModel`](@ref QuantumModels.QuantumModel) instead provides 
multiple potential energy surfaces with couplings between them. As implied by the name,
these are in the diabatic representation.
If the desired model does not fall under either of these branches, a new abstract type
should be created.

### Minor branches

Under the two main branches there also exists some specialised abstract types that are
useful in some cases, such as when using many electronic states, or when implementing
extra functions is required. See the docstrings for more info:

- [`ClassicalFrictionModel`](@ref FrictionModels.ClassicalFrictionModel)
- [`QuantumFrictionModel`](@ref FrictionModels.QuantumFrictionModel)

## The N-dimensional, M-atom convention

Every [`Model`](@ref NQCModels.Model) must work for an arbitrary number of atoms `M`
and degrees of freedom per atom `N = ndofs(model)`. Concretely, `R` is always an
`N * M` shaped `AbstractMatrix` (one column per atom), and:

- `potential(model, R)`/`potential!(model, V, R)` always return/fill a value sized for
  the *whole* system (a scalar for [`ClassicalModel`](@ref ClassicalModels.ClassicalModel),
  a single `nstates * nstates` `Hermitian` for [`QuantumModel`](@ref QuantumModels.QuantumModel)).
- `derivative(model, R)`/`derivative!(model, D, R)` always return/fill an `N * M` shaped
  array with one entry per (degree of freedom, atom) pair.
- `potential!`/`derivative!` **overwrite** `V`/`D` completely on every call rather than
  accumulating across calls — callers (e.g. `NQCCalculators.jl`) rely on this to reuse the
  same pre-allocated arrays every timestep without manually resetting them first.

Every model implements `potential`/`potential!`/`derivative!` directly against
`R::AbstractMatrix` — there is no per-atom "hook" method and no automatic fallback for
multi-atom systems. If a model does not explicitly define an interaction between different
atoms, its implementation should say so explicitly, typically by summing/broadcasting over
`eachcol(R)`/`axes(R, 2)` internally (e.g. [`Morse`](@ref)). Models that genuinely couple
different atoms (e.g. [`DiatomicHarmonic`](@ref), which couples exactly its first two atoms
and leaves any further atoms independent/non-interacting) or that have heterogeneous
per-atom parameters (e.g. [`Harmonic`](@ref)) are written the same way — see the
[`ClassicalModel`](@ref ClassicalModels.ClassicalModel) and
[`QuantumModel`](@ref QuantumModels.QuantumModel) docstrings for worked examples.

### Backward-compatible scalar entry

For models where `ndofs(model) == 1`, the position of a single atom can also be passed
directly as a plain number, e.g. `potential(model, 1.0)`, rather than as a `1×1` matrix.
This is handled generically by convenience methods defined on [`Model`](@ref NQCModels.Model)
that wrap the number in a `1×1` matrix and redirect to the matrix-based methods, so
individual models do not need to implement this themselves.

## Example implementations

To implement a new model, you first select the abstract type, where you should first take a look
at the docstring for the abstract type (click on the links above). There, a list of the functions 
that need to be implemented along with an example implementation are provided.

For further examples, you can also take a look into the source code of the `NQCModels.jl`
package to see how the analytic models have been implemented. 
If you have any issues or questions about implementing a new model, open up an issue on
Github and we can work together to resolve the problem.

## GPU Support Implementation

GPU support for models is provided through package extensions for CUDA.jl and Metal.jl. When implementing a new model, GPU support can be added following these guidelines:

### Extension Structure

GPU implementations are located in:
- `ext/NQCModelsCUDAExt.jl` - CUDA.jl (NVIDIA GPUs) support
- `ext/NQCModelsMetalExt.jl` - Metal.jl (Apple Silicon) support

### Type-Based Dispatch

GPU methods use Julia's multiple dispatch to activate when GPU arrays are passed:

```julia
# CPU version (in main package)
function NQCModels.potential(model::YourModel, R::AbstractMatrix)
    # ... CPU implementation
end

# GPU version (in extension)
function NQCModels.potential(model::YourModel, R::CUDA.CuArray)
    # ... GPU implementation using CUDA operations
end
```

### Implementation Guidelines

1. **Implement in-place methods only**: Only `potential!` and `derivative!` need GPU versions. The allocating versions (`potential`, `derivative`) should already be defined by model-unspecific methods unless a deviation from normal behaviour is needed.

2. **Match return types to input**: Return GPU arrays when given GPU input:
   ```julia
   function NQCModels.potential(model::YourModel, R::CUDA.CuArray)
       # Compute on GPU
       V_gpu = ...  # CUDA.CuArray
       return V_gpu
   end
   ```

3. **No Hermitian wrappers on GPU**: For quantum models, return plain GPU matrices rather than `Hermitian{T, CuArray}` to avoid complexity:
   ```julia
   # Return CuMatrix{Float64}, not Hermitian{Float64, CuMatrix{Float64}}
   return CUDA.CuArray([v11 v12; v21 v22])
   ```

4. **Add reference comments to CPU version**: Include line number references in the CPU implementation:
   ```julia
   # GPU extensions: ext/NQCModelsCUDAExt.jl:123, ext/NQCModelsMetalExt.jl:456
   function NQCModels.potential(model::YourModel, R::AbstractMatrix)
       # ... CPU implementation
   end
   ```

5. **Implement for both CUDA and Metal**: Add methods to both extension files to support all GPU platforms.

### Example GPU Implementation

Here's a complete example for a simple harmonic model:

```julia
# In ext/NQCModelsCUDAExt.jl
function NQCModels.potential(model::Harmonic, r::CUDA.CuArray)
    m = model.m
    ω = model.ω
    r₀ = model.r₀
    
    # Use CUDA operations
    return sum((m / 2) .* ω^2 .* (r .- r₀).^2)
end

function NQCModels.derivative(model::Harmonic, r::CUDA.CuArray)
    m = model.m
    ω = model.ω
    r₀ = model.r₀
    
    # Return CuArray
    return m .* ω^2 .* (r .- r₀)
end

# Repeat similar implementation in ext/NQCModelsMetalExt.jl
# using Metal.MtlArray instead of CUDA.CuArray
```

### Limitations

Not all models can support GPU arrays. Models that depend on:
- External libraries without GPU support (e.g., BSplineKit for `CSVModel_1D`)
- Python interop (e.g., JuLIP models)
- External calculators (e.g., AtomsCalculators models)

cannot provide GPU implementations.

### Testing GPU Support

Basic test pattern for GPU implementations:

```julia
using NQCModels
using CUDA
using Test

model = YourModel()
R_cpu = [1.0 2.0]
R_gpu = CUDA.CuArray(R_cpu)

V_cpu = potential(model, R_cpu)
V_gpu = potential(model, R_gpu)

@test V_cpu ≈ Array(V_gpu)  # Transfer GPU result back for comparison

D_cpu = derivative(model, R_cpu)
D_gpu = derivative(model, R_gpu)

@test D_cpu ≈ Array(D_gpu)
```

### Updating Documentation

When adding GPU support to a model:
1. Update the supported models list in `docs/src/NQCModels/gpu_support.md`
2. Add line number reference comments to the CPU implementation
3. Ensure both CUDA and Metal extensions are updated consistently

For more information on using GPU-accelerated models, see [GPU Acceleration](@ref). 
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
