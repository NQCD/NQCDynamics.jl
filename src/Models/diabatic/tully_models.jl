export TullyModelOne
export TullyModelTwo

"""
The 2-state Tully model 1
"""
struct TullyModelOne <: DiabaticModel
    n_states::UInt8
    A::Float64
    B::Float64
    C::Float64
    D::Float64
end

TullyModelOne(A=0.01, B=1.6, C=0.005, D=1.0) = TullyModelOne(2, A, B, C, D)

function potential!(model::TullyModelOne, V::Hermitian, R::AbstractMatrix)
    q = R[1]
    if q > 0
        V[1,1] = model.A * (1 - exp(-model.B*q))
    else
        V[1,1] = -model.A * (1 - exp(model.B*q))
    end
    V[2,2] = -V[1,1]
    V.data[1,2] = model.C * exp(-model.D*q^2) # sets both off-diagonals as V is Hermitian
end

function derivative!(model::TullyModelOne, derivative::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    q = R[1]
    derivative[1][1,1] = model.A * model.B * exp(-model.B * abs(q))
    derivative[1][2,2] = -derivative[1][1,1]
    derivative[1].data[1,2] = -2 * model.C * model.D * q * exp(-model.D*q^2)
end

struct TullyModelTwo <: DiabaticModel
    n_states::UInt8
    A::Float64
    B::Float64
    C::Float64
    D::Float64
    E::Float64
end

TullyModelTwo(A=0.1, B=0.28, C=0.015, D=0.06, E=0.05) = TullyModelTwo(2, A, B, C, D, E)

function potential!(model::TullyModelTwo, V::Hermitian, R::AbstractMatrix)
    q = R[1]
    V[1,1] = 0
    V[2,2] = -model.A*exp(-model.B*q^2) + model.E
    V.data[1,2]= model.C * exp(-model.D*q^2)
end

function derivative!(model::TullyModelTwo, derivative::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    q = R[1]
    derivative[1][1,1] = 0
    derivative[1][2,2] = 2*model.A*model.B*q*exp(-model.B*q^2)
    derivative[1].data[1,2] = -2*model.C*model.D*q*exp(-model.D*q^2)
end
