using NonadiabaticMolecularDynamics.Models
using Plots
using PeriodicTable
using UnitfulAtomic
using Unitful

N = 50
mass = austrip(elements[:H].atomic_mass)
D = austrip(2u"eV")
a = austrip(1u"Å")
B = sqrt(D / 4π)
β = D
model = Scattering1D(N=N, D=D, B=B, β=β, α=2β, a=a)

plot(-1.0:0.05:3, model)
