export Subotnik_A

@doc raw"""
    Model A according to Coffmann & Subotnik, Phys.Chem.Chem.Phys., 20, 9847 (2018).

    V_0 =  0.5*m*ω^2 r^2
    V_1 =  0.5*m*ω^2 (r-g)^2 + ΔG
    Molecule-metal coupling:
    Γ(x) = Γ_1(Γ_0 (K(r-d)^2)/(1+K(r-d)^2))
    d is the minimum of the PES, d = g/2 +ΔG/(mω^2g) 
    In the original model, K = 0.1 and Γ_0=0.01; Γ_1 is used as an energy unit
    The numbers given as default correspond to the one listed in the paper.
"""

struct Subotnik_A <: DiabaticModel
    n_states::UInt8
    m::Float64
    om::Float64
    g::Float64
    DG::Float64
    G1::Float64
    G0::Float64
    K:: Float64
    d::Float64
end

#function Subotnik_A(
#    ;N=2,
#    m=2000.0,
#    om=0.00002,
#    g=20.6097,
#   DG=-0.0038,
#    G1=0.0001,
#    G0=0.01,
#    K= 0.1,
#    d=8.0
#    )
Subotnik_A(m=2000.0, om=0.0002, g=20.6097, DG=-0.0038, G1=0.0001, G0=0.01, K=0.1, d=8.0) = Subotnik_A(2,m, om, g, DG, G1, G0, K, d)
#end

"""
Diabatic potenial matrix
"""
function potential!(model::Subotnik_A, V::Hermitian, R::AbstractMatrix)
    # Update d
    d = model.g/2 + model.DG/(model.m*model.om^2*model.g)


    # Higher potential well
    V_0(r) = (model.m*model.om^2*r^2)/2.0
    # Loweer energy potential well
    V_1(r) = (model.m*model.om^2*(r-model.g)^2)/2.0 + model.DG
    println(R[1]," " , V_0(R[1]), " ", V_1(R[1]))
    # weighted metal-molecule coupling between the diabates
    G(r)=model.G1*(model.G0 + (model.K*(r - model.d)^2)/(1.0 + model.K*(r -model.d)^2))
    
    
    V = Hermitian(zeros(2,2))
    # Set up the Hamiltonian
    V[1,1] = V_0(R[1])
    V[2,2] = V_1(R[1])
    V.data[1,2] = G(R[1])
    V.data[2,1] = G(R[1])
end

"""
Elementwise position derivative of the above potential
"""
function derivative!(model::Subotnik_A,D::AbstractMatrix{<:Hermitian}, R::AbstractArray)
    dV0dr(r) = model.m*model.om^2*r
    dV1dr(r) = model.m*model.om^2*(r-model.g)
    dGdr(r)  = model.G1((2*model.K*(r-model.d)*(1+K*(r-model.d)^2)-2*model.K^2*(r-model.d)^2)/
                         (1+model.K*(r-model.d)^2)^2)
    D[1][1,1] = dV0dr([R[1]])
    D[1][2,2] = dV1dr([R[1]])
    D[1][1,2] = dGdr([R[1]])
    D[1][2,1] = dGdr([R[1]])
end