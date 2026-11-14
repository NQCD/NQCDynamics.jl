
abstract type AnanthModel <: QuantumModel end

NQCModels.ndofs(::AnanthModel) = 1
NQCModels.nstates(::AnanthModel) = 2

"""
    AnanthModelOne(a=0.01, b=1.6, c=0.005, d=1.0)

Ananth's simple avoided crossing model (similar to Tully's first model) from [J. Chem. Phys. 127, 084114 (2007)](http://dx.doi.org/10.1063/1.2759932).
"""
Parameters.@with_kw struct AnanthModelOne{A,B,C,D} <: AnanthModel
    a::A = 0.01
    b::B = 1.6
    c::C = 0.005
    d::D = 1.0
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:204, ext/NQCModelsMetalExt.jl:204
function NQCModels.potential(model::AnanthModelOne, q::AbstractMatrix)
    Parameters.@unpack a, b, c, d = model
    V11 = sum(a .* tanh.(b .* q))
    V22 = sum(.-a .* tanh.(b .* q))
    V12 = sum(c .* exp.(.-d .* q .^ 2)) # sets both off-diagonals as V is Hermitian
    return Hermitian([V11 V12; V12 V22])
end

function NQCModels.potential!(model::AnanthModelOne, V::Hermitian, q::AbstractMatrix)
    Parameters.@unpack a, b, c, d = model
    V11 = sum(a .* tanh.(b .* q))
    V22 = sum(.-a .* tanh.(b .* q))
    V12 = sum(c .* exp.(.-d .* q .^ 2)) # sets both off-diagonals as V is Hermitian
    V.data .= [V11 V12; V12 V22]
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:216, ext/NQCModelsMetalExt.jl:216

function NQCModels.derivative!(model::AnanthModelOne, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    Parameters.@unpack a, b, c, d = model
    for dof in axes(R, 1)
        for atom in axes(R, 2)
            q = R[dof, atom]
            D11 = a * b * (1 - tanh(b*q)^2)
            D22 = -D11
            D12 = -2 * c * d * q * exp(-d*q^2)
            D[dof,atom].data .= [D11 D12; D12 D22]
        end
    end
end

"""
    AnanthModelTwo(a=0.04, b=0.01, c=0.005, d=1.0, e=0.7, f=1.6)

Ananth's asymmetric model from [J. Chem. Phys. 127, 084114 (2007)](http://dx.doi.org/10.1063/1.2759932).
"""
Parameters.@with_kw struct AnanthModelTwo{A,B,C,D,E,F} <: AnanthModel
    a::A = 0.04
    b::B = 0.01
    c::C = 0.005
    d::D = 1.0
    e::E = 0.7
    f::F = 1.6
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:229, ext/NQCModelsMetalExt.jl:229
function NQCModels.potential(model::AnanthModelTwo, Q::AbstractMatrix)
    Parameters.@unpack a, b, c, d, e, f = model
    V11 = a .* tanh.(f .* Q)
    V22 = -b .* tanh.(f .* Q)
    V12 = c .* exp.(-d .* (Q .+ e) .^ 2)
    return Hermitian([sum(V11) sum(V12); sum(V12) sum(V22)])
end

function NQCModels.potential!(model::AnanthModelTwo, V::Hermitian, Q::AbstractMatrix)
    Parameters.@unpack a, b, c, d, e, f = model
    V11 = a .* tanh.(f .* Q)
    V22 = -b .* tanh.(f .* Q)
    V12 = c .* exp.(-d .* (Q .+ e) .^ 2)
    V.data .= [sum(V11) sum(V12); sum(V12) sum(V22)]
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:241, ext/NQCModelsMetalExt.jl:241
function NQCModels.derivative!(model::AnanthModelTwo, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    Parameters.@unpack a, b, c, d, e, f = model
    for dof in axes(R, 1)
        for atom in axes(R, 2)
            q = R[dof, atom]
            D11 = a * f * (1 - tanh(f*q)^2)
            D22 = -b * f * (1 - tanh(f*q)^2)
            D12 = -2 * c * d * (q+e) * exp(-d*(q+e)^2)

            D[dof,atom].data .= [D11 D12; D12 D22]
        end
    end
end
