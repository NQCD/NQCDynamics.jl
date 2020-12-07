function calc_NJ!(enj,μ,AM,rmin,rmax,V_ref,n,J)
    using FastGaussQuadrature
    #    CALCULATE VIBRATIONAL AND ROTATIONAL QUANTUM NUMBERS FOR
    #    A PRODUCT DIATOM
 
    #     FIND THE ROTATIONAL QUANTUM NUMBER J FROM THE ROTATIONAL
    #     ANGULAR MOMENTUM AM IN UNITS OF ħ.
 
    n=0
    J=0.5*(-1+√(1+4*AM²))
    #    FIND THE VIBRATIONAL QUANTUM NUMBER n, USING AN INVERSION
    #    OF THE SEMICLASSICAL RYDBERG-KLEIN-REES (RKR) APPROACH.
    #
    #    INITIALIZE SOME VARIABLES FOR THE DIATOM.
 
    T3=Q(2,3)-Q(1,3)
    T2=Q(2,2)-Q(1,2)
    T1=Q(2,1)-Q(1,1)
    RO=√(T1*T1+T2*T2+T3*T3)
 
    H,T,V = diatomic.calc_energy(bondlength=R0)
    VEFF=V-V_ref+0.5*AM²*ħ²/(μ*RO²)
 
    # DETERMINE BOUNDARIES OF THE SEMICLASSICAL INTEGRAL
    while (VEFF<enj and RZ<50.0)
        RZ=RZ+0.001
        H,T,V = diatomic.calc_energy(bondlength=RZ)
        VEFF=V-V_ref+0.5*AM²*ħ²/(μ*RZ²)
    end
    rmax=RZ

    RZ=R0
    while (VEFF<enj)
       RZ=RZ-0.001
       H,T,V = diatomc.calc_energy(bondlength=RZ)
       VEFF=V-V_ref+0.5*AM²*ħ²/(μ*RZ²)
    end
    rmin=RZ
 
    # EVALUATE THE SEMICLASSICAL INTEGRAL BY GAUSSIAN QUADRATURE
 
    #Gauss-Legendre quadrature Parameters
    n_nodes = 50
    nodes, weights = gausslegendre(n_nodes)
    ASUM=0.0
    ΔR=(rmax-rmin)/50
    do i=1,n_nodes
       RZ=ΔR*(i-1)
       R0 = rmin+RZ
       H,T,V = diatomic.calc_energy(bondlength=R0)
       VEFF=(V-V_ref)+0.5*AM²*ħ²/(μ*RZ²)
       if (enj>VEFF)
          ASUM=ASUM+weights[i]*√(enj-VEFF)
       end
    end

    n=√(8.0*μ)*ASUM/2π/ħ
    n=n-0.5
 
 end