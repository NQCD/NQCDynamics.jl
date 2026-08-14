```@setup logging
@info "Expanding src/NQCModels/gpu_support.md..."
start_time = time()
```

# GPU Acceleration

GPU array support is available for many models in NQCModels.jl via package extensions for CUDA.jl (NVIDIA GPUs) and Metal.jl (Apple Silicon GPUs). This allows potential and derivative calculations to be performed directly on GPU arrays without explicit data transfers.

## Quick Start

### CUDA (NVIDIA GPUs)

```julia
using NQCModels
using CUDA

# Create a model
model = TullyModelOne()

# Use GPU arrays for positions
R_gpu = CUDA.CuArray([1.0])

# Evaluate on GPU - results stay on GPU
V = potential(model, R_gpu)  # Returns scalar
D = derivative(model, R_gpu)  # Returns CuArray

# Transfer back to CPU if needed
V_cpu = Array(V)
D_cpu = Array(D)
```

### Metal (Apple Silicon)

```julia
using NQCModels
using Metal

# Create a model
model = Harmonic(m=1.0, ω=1.0, r₀=0.0)

# Use GPU arrays for positions
R_gpu = Metal.MtlArray([1.0 2.0 3.0])

# Evaluate on GPU
V = potential(model, R_gpu)  # Returns scalar
D = derivative(model, R_gpu)  # Returns MtlArray
```

## Supported Models

### Classical Models

All basic classical models support GPU arrays:
- [`Free`](@ref)
- [`Harmonic`](@ref)
- [`Morse`](@ref)
- [`Logistic`](@ref)
- [`DiatomicHarmonic`](@ref)
- [`DarlingHollowayElbow`](@ref)

### Quantum Models

The following quantum models support GPU arrays:
- Tully models: [`TullyModelOne`](@ref), [`TullyModelTwo`](@ref), [`TullyModelThree`](@ref)
- Ananth models: [`AnanthModelOne`](@ref), [`AnanthModelTwo`](@ref)
- [`DoubleWell`](@ref)
- [`MiaoSubotnik`](@ref)
- [`GatesHollowayElbow`](@ref)
- [`ThreeStateMorse`](@ref)
- [`ErpenbeckThoss`](@ref)
- [`AndersonHaldane`](@ref)
- [`OuyangModelOne`](@ref)
- [`Scattering1D`](@ref)

### Composite Models

The following composite models automatically support GPU when their underlying models do:
- [`AveragedPotential`](@ref)
- [`AndersonHolstein`](@ref)
- [`WideBandBath`](@ref)
- [`AdiabaticStateSelector`](@ref)

## Implementation Details

### Type-Based Dispatch

GPU support works through Julia's multiple dispatch mechanism. When you pass a GPU array to `potential` or `derivative`, the GPU-specific implementation is automatically selected:

```julia
# CPU version (existing code)
R_cpu = [1.0]
V_cpu = potential(model, R_cpu)  # Uses CPU implementation

# GPU version (automatic when using GPU arrays)
R_gpu = CUDA.CuArray([1.0])
V_gpu = potential(model, R_gpu)  # Uses GPU implementation
```

### Return Types

Functions return arrays matching the input type:
- **CPU input** → CPU output
- **GPU input** → GPU output

For quantum models with GPU arrays, matrices are returned as plain GPU arrays rather than `Hermitian`-wrapped types to avoid complexity.

### In-Place vs Allocating

Currently, only the allocating versions (`potential`, `derivative`) have GPU support. The in-place versions (`potential!`, `derivative!`) remain CPU-only. This design:
- Keeps existing CPU code unchanged
- Avoids complications with GPU array mutation
- Relies on GPU compiler optimizations for performance

## Performance Tips

### When to Use GPU Arrays

GPU arrays are most beneficial for:
- **Batch calculations**: Evaluating many positions simultaneously
- **Large systems**: Models with many degrees of freedom or states
- **Iterative workflows**: When data stays on GPU across multiple evaluations

GPU arrays may be **slower** for:
- Single point evaluations (transfer overhead dominates)
- Small models (insufficient parallelism)
- Models with complex control flow

### Minimizing Transfers

Keep data on the GPU throughout your workflow:

```julia
using NQCModels
using CUDA

model = TullyModelOne()

# Bad: Transfer for each evaluation
for i in 1:100
    R_cpu = [rand()]
    R_gpu = CUDA.CuArray(R_cpu)  # Transfer to GPU
    V = potential(model, R_gpu)
    V_cpu = Array(V)  # Transfer back to CPU
end

# Good: Keep data on GPU
R_gpu = CUDA.CuArray([0.0])
results = CUDA.CuArray{Float64}(undef, 100)
for i in 1:100
    R_gpu[1] = rand()  # Update on GPU
    results[i] = potential(model, R_gpu)
end
results_cpu = Array(results)  # Single transfer at end
```

### Composite Models

Be aware that some composite models (e.g., `ErpenbeckThoss`, `AndersonHaldane`) may internally call CPU implementations of underlying models, which can cause unexpected transfers. Monitor performance when using these models.

## Limitations

The following models do **not** support GPU arrays:
- `CSVModel_1D`: Uses BSplineKit which doesn't support GPU arrays
- JuLIP models: Require Python interop
- AtomsCalculators models: External calculator dependencies

## Requirements

- **Julia**: 1.9 or later (for package extension support)
- **CUDA.jl**: v5.x (for NVIDIA GPUs)
- **Metal.jl**: v1.x (for Apple Silicon)
- Appropriate GPU drivers installed

GPU support is provided through package extensions and is automatically loaded when you use `using CUDA` or `using Metal` alongside `NQCModels`.

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
