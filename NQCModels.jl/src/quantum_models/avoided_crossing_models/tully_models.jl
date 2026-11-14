
abstract type TullyModel <: QuantumModel{Diabatic} end

NQCModels.ndofs(::TullyModel) = 1
NQCModels.nstates(::TullyModel) = 2

"""
    TullyModelOne(a=0.01, b=1.6, c=0.005, d=1.0)

Tully's simple avoided crossing model from [J. Chem. Phys. 93, 1061 (1990)](https://doi.org/10.1063/1.459170).
"""
Parameters.@with_kw struct TullyModelOne{A,B,C,D} <: TullyModel
    a::A = 0.01
    b::B = 1.6
    c::C = 0.005
    d::D = 1.0
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:105, ext/NQCModelsMetalExt.jl:105

function NQCModels.potential!(model::TullyModelOne, V::Hermitian, R::AbstractMatrix)
    Parameters.@unpack a, b, c, d = model
    V11 = sum(R) do q
        q > 0 ? a * (1 - exp(-b*q)) : -a * (1 - exp(b*q))
    end
    V12 = sum(q -> c * exp(-d*q^2), R) # sets both off-diagonals as V is Hermitian
    V.data .= [V11 V12; V12 -V11]
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:119, ext/NQCModelsMetalExt.jl:119

function NQCModels.derivative!(model::TullyModelOne, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    Parameters.@unpack a, b, c, d = model
    @inbounds for dof in axes(R, 1)
        @inbounds for atom in axes(R, 2)
            q = R[dof, atom]
            D11 = a * b * exp(-b * abs(q))
            D22 = -D11
            D12 = -2 * c * d * q * exp(-d*q^2)
            D[dof,atom].data .= [D11 D12; D12 D22]
        end
    end
end

"""
    TullyModelTwo(a=0.1, b=0.28, c=0.015, d=0.06, e=0.05)

Tully's dual avoided crossing model from [J. Chem. Phys. 93, 1061 (1990)](https://doi.org/10.1063/1.459170).
"""
Parameters.@with_kw struct TullyModelTwo{A,B,C,D,E} <: TullyModel
    a::A = 0.1
    b::B = 0.28
    c::C = 0.015
    d::D = 0.06
    e::E = 0.05
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:130, ext/NQCModelsMetalExt.jl:130

function NQCModels.potential!(model::TullyModelTwo, V::Hermitian, R::AbstractMatrix)
    Parameters.@unpack a, b, c, d, e = model
    V22 = sum(q -> -a*exp(-b*q^2) + e, R)
    V12 = sum(q -> c * exp(-d*q^2), R)
    V.data .= [0 V12; V12 V22]
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:140, ext/NQCModelsMetalExt.jl:140

function NQCModels.derivative!(model::TullyModelTwo, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    Parameters.@unpack a, b, c, d  = model
    @inbounds for dof in axes(R, 1)
        @inbounds for atom in axes(R, 2)
            q = R[dof, atom]
            D11 = 0
            D22 = 2*a*b*q*exp(-b*q^2)
            D12 = -2*c*d*q*exp(-d*q^2)

            D[dof,atom].data .= [D11 D12; D12 D22]
        end
    end
end

"""
    TullyModelThree(a=0.0006, b=0.1, c=0.9)

Tully's extended coupling with reflection model from [J. Chem. Phys. 93, 1061 (1990)](https://doi.org/10.1063/1.459170).
"""
Parameters.@with_kw struct TullyModelThree{A,B,C} <: TullyModel
    a::A = 0.0006
    b::B = 0.1
    c::C = 0.9
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:151, ext/NQCModelsMetalExt.jl:151
function NQCModels.potential!(model::TullyModelThree, V::Hermitian, R::AbstractMatrix)
    Parameters.@unpack a, b, c = model
    n = length(R)
    V12 = sum(R) do q
        q > 0 ? b * (2 - exp(-c*q)) : b * exp(c*q)
    end
    V.data .= [a*n V12; V12 -a*n]
end

function NQCModels.derivative!(model::TullyModelThree, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    Parameters.@unpack a, b, c = model
    @inbounds for dof in axes(R,1)
        @inbounds for atom in axes(R, 2)
            q = R[dof, atom]
            D11 = 0
            D22 = 0
            D12 = q > 0 ? b * c * exp(-c*q) : b * c * exp(c*q)
            D[dof,atom].data .= [D11 D12; D12 D22]
        end
    end
end
