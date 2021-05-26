export TullyModelOne
export TullyModelTwo

"""
    TullyModelOne(a=0.01, b=1.6, c=0.005, d=1.0)

Tully's simple avoided crossing model from [J. Chem. Phys. 93, 1061 (1990)](https://doi.org/10.1063/1.459170).
"""
@with_kw struct TullyModelOne{A,B,C,D} <: DiabaticModel
    n_states::UInt8 = 2
    a::A = 0.01
    b::B = 1.6
    c::C = 0.005
    d::D = 1.0
end

function potential!(model::TullyModelOne, V::Hermitian, R::AbstractMatrix)
    @unpack a, b, c, d = model
    q = R[1]
    if q > 0
        V[1,1] = a * (1 - exp(-b*q))
    else
        V[1,1] = -a * (1 - exp(b*q))
    end
    V[2,2] = -V[1,1]
    V.data[1,2] = c * exp(-d*q^2) # sets both off-diagonals as V is Hermitian
end

function derivative!(model::TullyModelOne, derivative::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    @unpack a, b, c, d = model
    q = R[1]
    derivative[1][1,1] = a * b * exp(-b * abs(q))
    derivative[1][2,2] = -derivative[1][1,1]
    derivative[1].data[1,2] = -2 * c * d * q * exp(-d*q^2)
end

"""
    TullyModelTwo(a=0.1, b=0.28, c=0.015, d=0.06, e=0.05)

Tully's dual avoided crossing model from [J. Chem. Phys. 93, 1061 (1990)](https://doi.org/10.1063/1.459170).
"""
@with_kw struct TullyModelTwo{A,B,C,D,E} <: DiabaticModel
    n_states::UInt8 = 2
    a::A = 0.1
    b::B = 0.28
    c::C = 0.015
    d::D = 0.06
    e::E = 0.05
end

function potential!(model::TullyModelTwo, V::Hermitian, R::AbstractMatrix)
    @unpack a, b, c, d, e = model
    q = R[1]
    V[1,1] = 0
    V[2,2] = -a*exp(-b*q^2) + e
    V.data[1,2]= c * exp(-d*q^2)
end

function derivative!(model::TullyModelTwo, derivative::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    @unpack a, b, c, d  = model
    q = R[1]
    derivative[1][1,1] = 0
    derivative[1][2,2] = 2*a*b*q*exp(-b*q^2)
    derivative[1].data[1,2] = -2*c*d*q*exp(-d*q^2)
end
