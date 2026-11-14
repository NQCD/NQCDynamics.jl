using Test
using NQCDynamics
using ComponentArrays
using Unitful

# A free ring polymer moving with constant centroid velocity through a periodic cell.
# All beads start at the same position with the same velocity, so the spring forces vanish
# and wrapping each bead individually is identical to wrapping the centroid.
atoms = Atoms([:H])
model = NQCModels.Free(3)
cell = PeriodicCell([10.0 0 0; 0 10.0 0; 0 0 10.0])
sim = RingPolymerSimulation{Classical}(atoms, model, 4; temperature = 100u"K", cell = cell)

r = fill(5.0, 3, 1, 4)
v = fill(0.1, 3, 1, 4)
z = ComponentVector(v = v, r = r)

solution = run_dynamics(sim, (0.0, 300.0), z; output = OutputPosition, dt = 0.1)
Rs = reduce(vcat, vec.(solution[:OutputPosition]))
@test maximum(Rs) > 10 # check the ring polymer leaves the cell

solution = run_dynamics(
    sim,
    (0.0, 300.0),
    z;
    output = OutputPosition,
    callback = DynamicsUtils.CellBoundaryCallback(),
    dt = 0.1,
)
Rs = reduce(vcat, vec.(solution[:OutputPosition]))
@test all(-1e-4 .<= Rs .<= 10 + 1e-4) # Check all beads remain inside
