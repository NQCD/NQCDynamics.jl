using LinearAlgebra

"""
    calculate_diatomic_energy!(sim::Simulation, bond_length::T; height::T=10.0,
                               normal_vector::Vector{T}=[0.0, 0.0, 1.0]) where {T}

Returns potential energy of diatomic with `bond_length` at `height` from surface.

Orients molecule parallel to the surface at the specified height, the surface is assumed
to intersect the origin.
This requires that the model implicitly provides the surface, or works fine without one.
"""
function calculate_diatomic_energy(sim::Simulation, bond_length::T;
                                   height::T=10.0,
                                   normal_vector::Vector{T}=[0.0, 0.0, 1.0]) where {T}
    #=
    This first line only needs to be calculated once for the system so it can be moved
    outside the function and passed as an argument.
    =#
    orthogonals = nullspace(reshape(normal_vector, 1, :)) # Plane parallel to surface
    R = normal_vector .* height .+ orthogonals .* bond_length/2

    Calculators.evaluate_potential!(sim.calculator, R)

    sim.calculator.potential[1]
end

function calculate_COM(R::Matrix, P::Matrix, atoms::Atoms)
    #calculate centre of mass (COM_R) from diatomic coordinates (R)
    #calculates centre of mass momenta (COM_R) from diatomic momenta (P)
    #COM_R = vector, COM_P = vector
    m = atoms.masses
    COM_R = (R[0,:]*m[0]+R[1,:]*m[1])/m[0]+m[1]
    COM_P = (P[0,:]*m[0]+P[1,:]*m[1])/m[0]+m[1]
    return COM_R, COM_P
end

function calculate_reduced_mass(atoms::Atoms)
    #calc reduced mass (μ)
    m = atoms.masses
    μ = m[1]*m[2]/(m[1]+m[2])
    μ
end

function calculate_bond_length(R::Matrix)
    #get bond length of diatomic
    T3=R[2,3]-R[1,3]
    T2=R[2,2]-R[1,2]
    T1=R[2,1]-R[1,1]
    RO=sqrt(T1^2+T2^2+T3^2)
    R0
end

function calculate_force_constant(sim::Simulation, bond_length::T;
                                height::T=10.0,
                                normal_vector::Vector{T}=[0.0, 0.0, 1.0]) where {T}
    using Statistics
    #Force constant (f) calcuulated from V(x) = 1/2 * f * x^2

    bond_lengths = collect(0.5:0.01:4.0)
    V = Vector(size(bond_lengths))

    for i in range(size(bond_lengths))
        V[i] = calculate_diatomic_energy(bond_length=bond_lengths[i])
    end

    f = Statistics.mean(V / bond_lengths^2) * 2
    f
    end
end

function calculate_quantum(R::Vector{T},P::Vector{T},Evib,AM,V_ref)
    using FastGaussQuadrature
    #    CALCULATE VIBRATIONAL AND ROTATIONAL QUANTUM NUMBERS FOR
    #    A PRODUCT DIATOM
 
    #     FIND THE ROTATIONAL QUANTUM NUMBER J FROM THE ROTATIONAL
    #     ANGULAR MOMENTUM AM IN UNITS OF ħ.
 
    n=0
    J=0.5*(-1+sqrt(1+4*AM^2))
    μ = calculate_reduced_mass()
    #    FIND THE VIBRATIONAL QUANTUM NUMBER n, USING AN INVERSION
    #    OF THE SEMICLASSICAL RYDBERG-KLEIN-REES (RKR) APPROACH.
    #
    #    INITIALIZE SOME VARIABLES FOR THE DIATOM.
 
    R_in = calculate_bond_length(R)
 
    H,T,V = calculate_diatomic_energy(bondlength=R_in)
    Veff=V-V_ref+0.5*AM^2*ħ^2/(μ*RO^2)
 
    # DETERMINE BOUNDARIES OF THE SEMICLASSICAL INTEGRAL
    while Veff<Evib & RZ<50.0
       RZ=RZ+0.001
       H,T,V = calculate_diatomic_energy(bondlength=RZ)
       Veff=V-V_ref+0.5*AM^2*ħ^2/(μ*RZ^2)
    end
    rmax=RZ
 
    RZ=R_in
    while Veff<Evib
       RZ=RZ-0.001
       H,T,V = calculate_diatomic_energy(bondlength=RZ)
       Veff=V-V_ref+0.5*AM^2*ħ^2/(μ*RZ^2)
    end
    rmin=RZ
    limts = [rmin,rmax]
    # EVALUATE THE SEMICLASSICAL INTEGRAL BY GAUSSIAN QUADRATURE
 
    #Gauss-Legendre quadrature Parameters
    n_nodes = 50
    nodes, weights = gausslegendre(n_nodes)
    ASUM=0.0
    ΔR=(rmax-rmin)/n_nodes
    for i=1:n_nodes
       R0 = rmin+(ΔR*(i-1))
       H,T,V = calculate_diatomic_energy(bondlength=R0)
       Veff=(V-V_ref)+0.5*AM^2*ħ^2/(μ*RZ^2)
       if (Evib>Veff)
          ASUM=ASUM+weights[i]*sqrt(Evib-Veff)
       end
    end
 
    n=sqrt(8.0*μ)*ASUM/2π/ħ
    n=n-0.5
 
    n,J,limits
end

#Set rotational / angular momenta from selected rotational state
#NEEDS WORK
function calculate_ROTN(AM,EROT)
    #Prior name = ROTN
 
    # C         CALCULATE ANGULAR MOMENTUM, MOMENT OF INERTIA TENSOR,
    # C         ANGULAR VELOCITY, AND ROTATIONAL ENERGY
 
    # C         CALCULATE ANGULAR MOMENTUM.  THE CENTER OF MASS COORDINATES
    # C         QQ AND MOMENTA PP COME FROM SUBROUTINE CENMAS THROUGH COMMON
    # C         BLOCK WASTE.

    #reminder: QQ , PP are pos, momenta for COM. W(N) are atom masses i think.
    #AM = angular momentum, C1 C2 etc are unit conversions.

    if (NAM=0) 
       do I=1,2
          J=L(I)
          J3=3*J
          J2=J3-1
          J1=J2-1
          AM(1)=AM(1)+(QQ(J2)*PP(J3)-QQ(J3)*PP(J2))
          AM(2)=AM(2)+(QQ(J3)*PP(J1)-QQ(J1)*PP(J3))
          AM(3)=AM(3)+(QQ(J1)*PP(J2)-QQ(J2)*PP(J1))
       end
       AM(4)=sqrt(AM(1)**2+AM(2)**2+AM(3)**2)
       
       AIXX=0.0
       do I=1,2
          J=3*L(I)+1
          SR=0.0
          do K=1,3
             SR=SR+QQ(J-K)**2
          end
          AIXX=AIXX+SR*W(L(I))
       end
       EROT=AM(4)^2/AIXX/2
       return
    end
    #         CALCULATE THE MOMENT OF INERTIA TENSOR

    do I=1,2
       J=L(I)
       J3=3*J
       J2=J3-1
       J1=J2-1
       AIXX=AIXX+W(J)*(QQ(J2)**2+QQ(J3)**2)
       AIYY=AIYY+W(J)*(QQ(J1)**2+QQ(J3)**2)
       AIZZ=AIZZ+W(J)*(QQ(J1)**2+QQ(J2)**2)
       AIXY=AIXY+W(J)*QQ(J1)*QQ(J2)
       AIXZ=AIXZ+W(J)*QQ(J1)*QQ(J3)
       AIYZ=AIYZ+W(J)*QQ(J2)*QQ(J3)
    end
    DET=AIXX*(AIYY*AIZZ-AIYZ*AIYZ)-AIXY*(AIXY*AIZZ+AIYZ*AIXZ)-
    *AIXZ*(AIXY*AIYZ+AIYY*AIXZ)
    #
    #         CALCULATE INVERSE OF THE INERTIA TENSOR
    #
    if (ABS(DET)>=0.01)
       UXX=(AIYY*AIZZ-AIYZ*AIYZ)/DET
       UXY=(AIXY*AIZZ+AIXZ*AIYZ)/DET
       UXZ=(AIXY*AIYZ+AIXZ*AIYY)/DET
       UYY=(AIXX*AIZZ-AIXZ*AIXZ)/DET
       UYZ=(AIXX*AIYZ+AIXZ*AIXY)/DET
       UZZ=(AIXX*AIYY-AIXY*AIXY)/DET
       #     CALCULATE ANGULAR VELOCITIES
       WX=UXX*AM(1)+UXY*AM(2)+UXZ*AM(3)
       WY=UXY*AM(1)+UYY*AM(2)+UYZ*AM(3)
       WZ=UXZ*AM(1)+UYZ*AM(2)+UZZ*AM(3)
    else
       #    CALCULATE ROTATIONAL ENERGY
       AIXX=0.0
       do I=1,2
          J=3*L(I)+1
          SR=0.0
          do K=1,3
                SR=SR+QQ(J-K)^2
          end
          AIXX=AIXX+SR*W(L(I))
       end
       EROT=AM(4)^2/AIXX/2
       RETURN
    end
    EROT=(WX*AM(1)+WY*AM(2)+WZ*AM(3))/2
    RETURN
end