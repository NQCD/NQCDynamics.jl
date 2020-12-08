function calc_nJ(R::Vector{T},P::Vector{T},enj,AM,V_ref)
    using FastGaussQuadrature
    #    CALCULATE VIBRATIONAL AND ROTATIONAL QUANTUM NUMBERS FOR
    #    A PRODUCT DIATOM
 
    #     FIND THE ROTATIONAL QUANTUM NUMBER J FROM THE ROTATIONAL
    #     ANGULAR MOMENTUM AM IN UNITS OF ħ.
 
    n=0
    J=0.5*(-1+√(1+4*AM^2))
    μ = diatomic_reduced_mass()
    #    FIND THE VIBRATIONAL QUANTUM NUMBER n, USING AN INVERSION
    #    OF THE SEMICLASSICAL RYDBERG-KLEIN-REES (RKR) APPROACH.
    #
    #    INITIALIZE SOME VARIABLES FOR THE DIATOM.
 
    #get bond length of diatomic
    T3=R[2,3]-R[1,3]
    T2=R[2,2]-R[1,2]
    T1=R[2,1]-R[1,1]
    RO=√(T1^2+T2^2+T3^2)
 
    H,T,V = diatomic.calc_energy(bondlength=R0)
    VEFF=V-V_ref+0.5*AM^2*ħ^2/(μ*RO^2)
 
    # DETERMINE BOUNDARIES OF THE SEMICLASSICAL INTEGRAL
    while VEFF<enj & RZ<50.0
        RZ=RZ+0.001
        H,T,V = diatomic.calc_energy(bondlength=RZ)
        VEFF=V-V_ref+0.5*AM^2*ħ^2/(μ*RZ^2)
    end
    rmax=RZ
 
    RZ=R0
    while VEFF<enj
       RZ=RZ-0.001
       H,T,V = diatomc.calc_energy(bondlength=RZ)
       VEFF=V-V_ref+0.5*AM^2*ħ^2/(μ*RZ^2)
    end
    rmin=RZ
    limts = [rmin,rmax]
    # EVALUATE THE SEMICLASSICAL INTEGRAL BY GAUSSIAN QUADRATURE
 
    #Gauss-Legendre quadrature Parameters
    n_nodes = 50
    nodes, weights = gausslegendre(n_nodes)
    ASUM=0.0
    ΔR=(rmax-rmin)/50
    for i=1:n_nodes
       R0 = rmin+(ΔR*(i-1))
       H,T,V = diatomic.calc_energy(bondlength=R0)
       VEFF=(V-V_ref)+0.5*AM^2*ħ^2/(μ*RZ^2)
       if (enj>VEFF)
          ASUM=ASUM+weights[i]*√(enj-VEFF)
       end
    end
 
    n=√(8.0*μ)*ASUM/2π/ħ
    n=n-0.5

    return n,J,limits
 end