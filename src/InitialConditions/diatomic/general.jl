using LinearAlgebra

function calculate_centre_of_mass(R::Matrix, atoms::Atoms)
    #calculate centre of mass from diatomic coordinates


end

function calculate_centre_of_mass_momenta(P::Matrix, atoms::Atoms)

end

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

function diatomic_reduced_mass()
    #calc reduced mass (μ)
    

    return μ
end

#generate large seed for all random number generation. Should be able to rerun with seed

function gen_seed()

    return seed
end