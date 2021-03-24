using NonadiabaticMolecularDynamics
using Plots
using Unitful
using UnitfulAtomic

atoms = Atoms([:H])
model = Models.FreeConstantFriction(1) # Friction in atomic units of mass / time
sim = Simulation{MDEF}(atoms, model; temperature=2000u"K")

"""
Temperature as a function of time

Must be defined to take atomic units and return atomic units
"""
function temperature_function(t)
    t_ps = au_to_ps(t)
    k = 1e-1
    T = exp(-k*t_ps) * 2000u"K"
    austrip(T)
end

# Uncomment this to use TwoTemperatureMDEF instead
# sim = Simulation{TwoTemperatureMDEF}(atoms, model, temperature_function)

# Initial velocity and positions
velocity = zeros(3,1)
position = zeros(3,1)
z = ClassicalDynamicals(velocity, position)

solution = Dynamics.run_trajectory(z, (0.0, 1e2u"ps"), sim, dt=1u"fs", output=(:energy, :position, :velocity))
plot(solution, :energy)
