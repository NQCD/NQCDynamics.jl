using NonadiabaticMolecularDynamics
using Plots
using PeriodicTable
using UnitfulAtomic
using Unitful

gr(framestyle=:box, palette=:mk_8, markersize=8, linewidth=5, markerstrokewidth=3, gridlinewidth=3)

N = 100
mass = austrip(elements[:H].atomic_mass)
D = austrip(2u"eV")
a = austrip(1u"Å")
B = sqrt(D / 4π)*1
β = D

# plot(-1.0:0.05:10, model)

rs = -0:0.01:2.5

plt = plot(legend=:topleft)
labels = [0.1, 1, 10, 50]
colors = [1, 5, 3, 4]
Bs = B .* labels
σs = austrip.([1e-4, 1e-3, 1e-2, 1e-1, 0.6, 1] .* u"eV")
Ns = [125, 150, 200]
for (i, B) in enumerate(Bs)
    model = Models.Scattering1D(N=N, D=D, B=B, β=β, α=2β, a=a)#, σ=σ)
    calculator = NonadiabaticMolecularDynamics.Calculators.Calculator(model, 1, 1)
    couplings = []
    @time for r in rs
        r = hcat(r)
        Calculators.update_electronics!(calculator, r)
        Calculators.evaluate_friction!(calculator, r)
        f = calculator.friction[1]
        f = f / mass
        f = au_to_ps_inv(f)
        push!(couplings, f)
    end
    plot!(rs, couplings, label="B = $(labels[i])x", color=colors[i])
end
xlims!(0, 2.5)
xlabel!("r /a₀")
ylabel!("γ(r) /ps⁻¹")

plt
