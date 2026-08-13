
"""
    ThreeStateMorse()

Three state morse potential referred to as Model IA here:
[J. Chem. Phys. 150, 244102 (2019)](https://doi.org/10.1063/1.5096276) 

Models IB and IC retain the same functional form and need only a change of parameters.
"""
Parameters.@with_kw struct ThreeStateMorse <: QuantumModel
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

NQCModels.ndofs(::ThreeStateMorse) = 1
NQCModels.nstates(::ThreeStateMorse) = 3

TSM_V_ii(x, d, α, r, c) = d * (1 - exp(-α*(x-r)))^2 + c
TSM_V_ij(x, a, α, r) = a * exp(-α*(x-r)^2)
function TSM_D_ii(x, d, α, r)
        ex = exp(-α*(x-r))
        return 2 * d * α * (ex - ex^2)
end
TSM_D_ij(x, a, α, r) = -2 * a * α * (x-r) * exp(-α*(x-r)^2)

function NQCModels.potential(model::ThreeStateMorse, R::AbstractMatrix)
    r = R[1]

    V11 = TSM_V_ii(r, model.d1, model.α1, model.r1, model.c1)
    V22 = TSM_V_ii(r, model.d2, model.α2, model.r2, model.c2)
    V33 = TSM_V_ii(r, model.d3, model.α3, model.r3, model.c3)

    V12 = TSM_V_ij(r, model.a12, model.α12, model.r12)
    V13 = TSM_V_ij(r, model.a13, model.α13, model.r13)
    V23 = TSM_V_ij(r, model.a23, model.α23, model.r23)

    return Hermitian([V11 V12 V13; V12 V22 V23; V13 V23 V33])
end

function NQCModels.potential!(model::ThreeStateMorse, V::Hermitian, R::AbstractMatrix)
    r = R[1]

    V11 = TSM_V_ii(r, model.d1, model.α1, model.r1, model.c1)
    V22 = TSM_V_ii(r, model.d2, model.α2, model.r2, model.c2)
    V33 = TSM_V_ii(r, model.d3, model.α3, model.r3, model.c3)

    V12 = TSM_V_ij(r, model.a12, model.α12, model.r12)
    V13 = TSM_V_ij(r, model.a13, model.α13, model.r13)
    V23 = TSM_V_ij(r, model.a23, model.α23, model.r23)

    V.data .= Hermitian([V11 V12 V13; V12 V22 V23; V13 V23 V33])
end

function NQCModels.derivative(model::ThreeStateMorse, R::AbstractMatrix)
    r = R[1]
    
    D11 = TSM_D_ii(r, model.d1, model.α1, model.r1)
    D22 = TSM_D_ii(r, model.d2, model.α2, model.r2)
    D33 = TSM_D_ii(r, model.d3, model.α3, model.r3)

    D12 = TSM_D_ij(r, model.a12, model.α12, model.r12)
    D13 = TSM_D_ij(r, model.a13, model.α13, model.r13)
    D23 = TSM_D_ij(r, model.a23, model.α23, model.r23)

    return Hermitian([D11 D12 D13; D12 D22 D23; D13 D23 D33])
end

function NQCModels.derivative!(model::ThreeStateMorse, D::Hermitian, R::AbstractMatrix)
    r = R[1]

    D11 = TSM_D_ii(r, model.d1, model.α1, model.r1)
    D22 = TSM_D_ii(r, model.d2, model.α2, model.r2)
    D33 = TSM_D_ii(r, model.d3, model.α3, model.r3)

    D12 = TSM_D_ij(r, model.a12, model.α12, model.r12)
    D13 = TSM_D_ij(r, model.a13, model.α13, model.r13)
    D23 = TSM_D_ij(r, model.a23, model.α23, model.r23)

    D.data .= Hermitian([D11 D12 D13; D12 D22 D23; D13 D23 D33])
end