using Interpolations

function DynamicsMethods.frozennuclei!(du, u, p::Tuple{Real,Simulation}, t)
    r = p[1]
    sim = p[2]
    Calcultors.update_electronics!(sim.calculator, r)
    DynamicsUtils.set_quantum_derivative!(du, u, sim) # 
end

function DynamicsMethods.frozennuclei!(du, u, p::Tuple{Interpolations.AbstractInterpolation,Simulation}, t)
    r = p[1](t)
    sim = p[2]
    Calcultors.update_electronics!(sim.calculator, r)
    DynamicsUtils.set_quantum_derivative!(du, u, sim)
end
