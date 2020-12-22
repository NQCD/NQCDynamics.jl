using Test
using LinearAlgebra
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.Calculators

cell = PeriodicCell(Matrix(I,3,3))
calc = PdH([:Pd, :H], cell, 10.0)

R = rand(3, 2)
@test evaluate_potential(calc, R) isa AbstractFloat
@test evaluate_derivative(calc, R) isa Matrix

calc = Free()
@test evaluate_potential(calc, R) == 0.0
@test evaluate_derivative(calc, R) == zeros(3, 2)

calc = Harmonic()
@test evaluate_potential(calc, R) == sum(R.^2)/2
@test evaluate_derivative(calc, R) == R