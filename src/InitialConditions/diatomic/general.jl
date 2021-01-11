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

function calc_COM(R::Matrix, P::Matrix, atoms::Atoms)
    #calculate centre of mass (COM_R) from diatomic coordinates (R)
    #calculates centre of mass momenta (COM_R) from diatomic momenta (P)
    #COM_R = vector, COM_P = vector
    m = atoms.masses
    COM_R = (R[0,:]*m[0]+R[1,:]*m[1])/m[0]+m[1]
    COM_P = (P[0,:]*m[0]+P[1,:]*m[1])/m[0]+m[1]
    return COM_R, COM_P
end

function calc_μ(atoms::Atoms)
    #calc reduced mass (μ)
    m = atoms.masses
    μ = m[1]*m[2]/(m[1]+m[2])
    return μ
end

function calc_force_constant(R::Matrix,atoms::Atoms)

    return f
end
# #generate large seed for all random number generation. Should be able to rerun with seed

# function gen_seed()

#     return seed
# end