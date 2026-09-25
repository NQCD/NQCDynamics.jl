using Test
using NQCDynamics
using Unitful, UnitfulAtomic
using Interpolations: linear_interpolation

@testset "Temperature - constant scalar" begin
    # Plain dimensionless number (assumed atomic units)
    T = Temperature(1.0e-3)
    @test T isa Temperature
    @test T(0.0) isa Real
    @test T(0.0) == 1.0e-3
    @test T(100.0) == 1.0e-3  # constant - independent of time

    # Unitful quantity
    T_k = Temperature(300u"K")
    @test T_k isa Temperature
    @test T_k(0.0) isa Real
    @test T_k(0.0) ≈ austrip(300u"K")
    @test T_k(50.0) ≈ austrip(300u"K")  # constant

    # Zero temperature
    T_zero = Temperature(0u"K")
    @test T_zero(0.0) == 0.0
end

@testset "Temperature - per-atom vector" begin
    # Plain dimensionless vector
    vals = [1.0e-3, 2.0e-3]
    T = Temperature(vals)
    @test T isa Temperature
    @test T(0.0) isa AbstractVector
    @test T(0.0) == vals
    @test T(99.0) == vals  # constant

    # Unitful vector
    T_k = Temperature([100u"K", 200u"K"])
    @test T_k(0.0) isa AbstractVector
    @test length(T_k(0.0)) == 2
    @test T_k(0.0)[1] ≈ austrip(100u"K")
    @test T_k(0.0)[2] ≈ austrip(200u"K")
end

@testset "Temperature - time-dependent function" begin
    # Function returning dimensionless number
    T = Temperature(t -> 1.0e-3 + t * 1.0e-6)
    @test T isa Temperature
    @test T(0.0) ≈ 1.0e-3
    @test T(10.0) ≈ 1.0e-3 + 10.0 * 1.0e-6

    # Function returning Unitful quantity
    T_k = Temperature(t -> 300u"K")
    @test T_k(0.0) isa Real
    @test T_k(0.0) ≈ austrip(300u"K")
end

@testset "Temperature - spline (times, values)" begin
    # Build from plain-number arrays (atomic units)
    times = [0.0, 10.0, 20.0]
    values = [1.0e-3, 2.0e-3, 3.0e-3]
    T = Temperature(times, values)
    @test T isa Temperature
    @test T(0.0) ≈ 1.0e-3
    @test T(10.0) ≈ 2.0e-3
    @test T(20.0) ≈ 3.0e-3
    # Interpolated value between knots
    @test T(5.0) ≈ 1.5e-3

    # Build from Unitful time arrays
    times_u = [0.0u"fs", 1.0u"fs", 2.0u"fs"]
    values_u = [100u"K", 200u"K", 300u"K"]
    T_u = Temperature(times_u, values_u)
    @test T_u isa Temperature
    @test T_u(austrip(0.0u"fs")) ≈ austrip(100u"K")
    @test T_u(austrip(1.0u"fs")) ≈ austrip(200u"K")
    @test T_u(austrip(2.0u"fs")) ≈ austrip(300u"K")
end

@testset "Temperature - get_temperature integration" begin
    # Test get_temperature method with Temperature objects directly
    T_obj = Temperature(100u"K")
    @test NQCDynamics.get_temperature(T_obj) isa Real
    @test NQCDynamics.get_temperature(T_obj) ≈ austrip(100u"K")
    @test NQCDynamics.get_temperature(T_obj, 5.0) ≈ austrip(100u"K")

    # Function-based Temperature object
    T_fn = Temperature(t -> 50u"K")
    @test NQCDynamics.get_temperature(T_fn) isa Real
    @test NQCDynamics.get_temperature(T_fn) ≈ austrip(50u"K")

    # Spline-based Temperature object
    times = austrip.([0.0u"fs", 1.0u"fs"])
    values = austrip.([100u"K", 200u"K"])
    T_spline = Temperature(times, values)
    @test NQCDynamics.get_temperature(T_spline, austrip(0.0u"fs")) ≈ austrip(100u"K")
    @test NQCDynamics.get_temperature(T_spline, austrip(1.0u"fs")) ≈ austrip(200u"K")
end
