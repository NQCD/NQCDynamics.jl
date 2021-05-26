
export ThreeStateMorse

"""
    ThreeStateMorse()

Three state morse potential referred to as Model IA here:
[J. Chem. Phys. 150, 244102 (2019)](https://doi.org/10.1063/1.5096276) 

Models IB and IC retain the same functional form and need only a change of parameters.
"""
@with_kw struct ThreeStateMorse <: DiabaticModel
    n_states::UInt8 = 3

    d1::Float64 = 0.02
    d2::Float64 = 0.02
    d3::Float64 = 0.003

    α1::Float64 = 0.4
    α2::Float64 = 0.65
    α3::Float64 = 0.65

    r1::Float64 = 4.0
    r2::Float64 = 4.5
    r3::Float64 = 6.0

    c1::Float64 = 0.02
    c2::Float64 = 0.0
    c3::Float64 = 0.02

    a12::Float64 = 0.005
    a13::Float64 = 0.005
    a23::Float64 = 0.0

    α12::Float64 = 32.0
    α13::Float64 = 32.0
    α23::Float64 = 0.0

    r12::Float64 = 3.40
    r13::Float64 = 4.97
    r23::Float64 = 0.0
end

function potential!(model::ThreeStateMorse, V::Hermitian, R::AbstractMatrix)
    V_ii(x, d, α, r, c) = d * (1 - exp(-α*(x-r)))^2 + c
    V_ij(x, a, α, r) = a * exp(-α*(x-r)^2)

    V[1,1] = V_ii(R[1], model.d1, model.α1, model.r1, model.c1)
    V[2,2] = V_ii(R[1], model.d2, model.α2, model.r2, model.c2)
    V[3,3] = V_ii(R[1], model.d3, model.α3, model.r3, model.c3)

    V.data[1,2] = V_ij(R[1], model.a12, model.α12, model.r12)
    V.data[1,3] = V_ij(R[1], model.a13, model.α13, model.r13)
    V.data[2,3] = V_ij(R[1], model.a23, model.α23, model.r23)
end

function derivative!(model::ThreeStateMorse, derivative::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)

    function D_ii(x, d, α, r)
        ex = exp(-α*(x-r))
        return 2 * d * α * (ex - ex^2)
    end

    D_ij(x, a, α, r) = -2 * a * α * (x-r) * exp(-α*(x-r)^2)

    derivative[1][1,1] = D_ii(R[1], model.d1, model.α1, model.r1)
    derivative[1][2,2] = D_ii(R[1], model.d2, model.α2, model.r2)
    derivative[1][3,3] = D_ii(R[1], model.d3, model.α3, model.r3)

    derivative[1].data[1,2] = D_ij(R[1], model.a12, model.α12, model.r12)
    derivative[1].data[1,3] = D_ij(R[1], model.a13, model.α13, model.r13)
    derivative[1].data[2,3] = D_ij(R[1], model.a23, model.α23, model.r23)
end
