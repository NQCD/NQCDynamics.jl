export MetallicChain
using LinearAlgebra
using BlockArrays

"""
The model supposed to emulate an adatom on a chain of metallic atoms.
Used in the paper below.
J. Phys. Chem. Lett. 2017, 8, 440−444
"""
struct MetallicChain <: AnalyticModel

    @add_standard_fields

    function MetallicChain(D=5.0::Float64, mass=2000.0::Float64, Ω=0.0028::Float64,
                           n=4::Integer, m=4::Integer,
                           d=1e-4::Float64, α=0.004::Float64,
                           β=0.03::Float64, γ=0.02::Float64, δ=0.01::Float64)

        alphas = fill(α, m)
        m_range = collect(0:m)
        betas = Diagonal(fill(β, m+1))
        gammas = Diagonal(fill(γ, m+1))
        deltas = Diagonal(fill(δ, m+1))

        V(R::Float64) = 0.5 * mass * Ω^2 * R^2
        Vkd(R::Float64, k::Int) = V(R - k*D)
        VkDdα(R::Float64, k::Int) = Tridiagonal(alphas, Vkd.(R .+ m_range.*d, k), alphas)

        dVdR(R::Float64) = mass * Ω * R
        dVdRkd(R::Float64, k::Int) = dVdR(R - k*D)
        dVdRkDdα(R::Float64, k::Int) = dVdRkd.(R .+ m_range.*d, k)

        function potential(q)::Hermitian
            b = BlockArray{Float64}(undef, fill(m+1, 2n+1), fill(m+1, 2n+1))
            for i=1:2n+1
                for j=1:2n+1
                    if i==j
                        setblock!(b, VkDdα(q, i-n-1), i, j)
                    elseif abs(i-j) == 1
                        setblock!(b, betas, i, j)
                    elseif abs(i-j) == 2
                        setblock!(b, gammas, i, j)
                    else
                        setblock!(b, deltas, i, j)
                    end
                end
            end

            return Hermitian(Array(b))
        end

        Deriv(q) = Hermitian(Diagonal(vcat([dVdRkDdα(q, i-n-1) for i=1:2n+1]...)))

        new(zero, zero, potential, Deriv)
    end
end
