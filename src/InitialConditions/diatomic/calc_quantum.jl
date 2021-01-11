function calc_nJ(R::Vector{T},P::Vector{T},Evib,AM,V_ref)
   using FastGaussQuadrature
   #    CALCULATE VIBRATIONAL AND ROTATIONAL QUANTUM NUMBERS FOR
   #    A PRODUCT DIATOM

   #     FIND THE ROTATIONAL QUANTUM NUMBER J FROM THE ROTATIONAL
   #     ANGULAR MOMENTUM AM IN UNITS OF ħ.

   n=0
   J=0.5*(-1+sqrt(1+4*AM^2))
   μ = calc_μ()
   #    FIND THE VIBRATIONAL QUANTUM NUMBER n, USING AN INVERSION
   #    OF THE SEMICLASSICAL RYDBERG-KLEIN-REES (RKR) APPROACH.
   #
   #    INITIALIZE SOME VARIABLES FOR THE DIATOM.

   R_in = calc_bond_length(R)

   H,T,V = calculate_diatomic_energy(bondlength=R_in)
   VEFF=V-V_ref+0.5*AM^2*ħ^2/(μ*RO^2)

   # DETERMINE BOUNDARIES OF THE SEMICLASSICAL INTEGRAL
   while VEFF<Evib & RZ<50.0
      RZ=RZ+0.001
      H,T,V = calculate_diatomic_energy(bondlength=RZ)
      VEFF=V-V_ref+0.5*AM^2*ħ^2/(μ*RZ^2)
   end
   rmax=RZ

   RZ=R_in
   while VEFF<Evib
      RZ=RZ-0.001
      H,T,V = calculate_diatomic_energy(bondlength=RZ)
      VEFF=V-V_ref+0.5*AM^2*ħ^2/(μ*RZ^2)
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
      VEFF=(V-V_ref)+0.5*AM^2*ħ^2/(μ*RZ^2)
      if (Evib>VEFF)
         ASUM=ASUM+weights[i]*sqrt(Evib-VEFF)
      end
   end

   n=sqrt(8.0*μ)*ASUM/2π/ħ
   n=n-0.5

   n,J,limits
end