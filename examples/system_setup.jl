using MDPrePostProcessing.Systems
using Unitful

x = [10, 0, 0]
y = [0, 10, 0]
z = [0, 0, 10]
cell = Cell([x y z], u"Å")  # Simulation cell

atom_types = [fill(:Pt, 10); :O; :C]  # Vector of all elements in the system
p = SystemParameters(cell, atom_types)  # Static system parameters

@show p.cell
@show p.atom_types
@show p.masses
@show p.n_atoms
@show p.n_beads

# Randomly select positions and momenta
R = rand(p.n_atoms)
P = rand(p.n_atoms)
z = Phasespace(R, P)

@show get_positions(z)
@show get_momenta(z)

system = System(p, z)