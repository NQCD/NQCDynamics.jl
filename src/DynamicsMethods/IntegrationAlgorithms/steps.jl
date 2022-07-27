using NQCDynamics: get_ring_polymer_temperature

function step_B!(v2, v1, dt, k)
    @.. broadcast=false v2 = muladd(dt, k, v1)
end
function step_A!(r2, r1, dt, v)
    @.. broadcast=false r2 = muladd(dt, v, r1)
end

function step_C!(v::AbstractArray{T,3}, r::AbstractArray{T,3}, cayley::Vector{<:AbstractMatrix}) where {T}
    for I in CartesianIndices(v)
        rtmp = cayley[I[3]][1,1] * r[I] + cayley[I[3]][1,2] * v[I]
        vtmp = cayley[I[3]][2,1] * r[I] + cayley[I[3]][2,2] * v[I]
        r[I] = rtmp
        v[I] = vtmp
    end
end

function step_C!(v::RingPolymerArray, r::RingPolymerArray, cayley::Vector{<:Matrix})
    for i in axes(r, 3)
        for j in RingPolymerArrays.quantumindices(v)
            for k in axes(r, 1)
                rtmp = cayley[i][1,1] * r[k,j,i] + cayley[i][1,2] * v[k,j,i]
                vtmp = cayley[i][2,1] * r[k,j,i] + cayley[i][2,2] * v[k,j,i]
                r[k,j,i] = rtmp
                v[k,j,i] = vtmp
            end
        end
    end

    for j in RingPolymerArrays.classicalindices(v)
        for k in axes(r, 1)
            rtmp = cayley[1][1,1] * r[k,j,1] + cayley[1][1,2] * v[k,j,1]
            vtmp = cayley[1][2,1] * r[k,j,1] + cayley[1][2,2] * v[k,j,1]
            r[k,j,1] = rtmp
            v[k,j,1] = vtmp
        end
    end
end

function step_O!(cache, integrator)
    @unpack t, dt, W, p, sqdt = integrator
    @unpack dutmp, flatdutmp, tmp1, tmp2, gtmp, noise, half, c1, c2 = cache

    Λ = gtmp
    σ = sqrt(get_temperature(p, t+dt*half)) ./ sqrt.(repeat(p.atoms.masses; inner=ndofs(p)))

    @.. noise = σ*W.dW[:] / sqdt

    γ, c = LAPACK.syev!('V', 'U', Λ) # symmetric eigen
    clamp!(γ, 0, Inf)
    for (j,i) in enumerate(diagind(c1))
        c1[i] = exp(-γ[j]*dt)
        c2[i] = sqrt(1 - c1[i]^2)
    end

    # equivalent to: flatdutmp .= c*c1*c'dutmp[:] + c*c2*c'noise
    copyto!(flatdutmp, dutmp)
    mul!(tmp1, transpose(c), flatdutmp)
    mul!(tmp2, c1, tmp1)
    mul!(flatdutmp, c, tmp2)
    
    mul!(tmp1, transpose(c), noise)
    mul!(tmp2, c2, tmp1)
    mul!(tmp1, c, tmp2)

    @.. flatdutmp += tmp1
    copyto!(dutmp, flatdutmp)
end

function step_O!(friction::MDEFCache, integrator, v, r, t)
    @unpack W, p, dt, sqdt = integrator
    @unpack flatvtmp, tmp1, tmp2, gtmp, noise, c1, c2, sqrtmass, σ = friction

    RingPolymerArrays.transform_from_normal_modes!(r, p.beads.transformation)
    RingPolymerArrays.transform_from_normal_modes!(v, p.beads.transformation)

    integrator.g(gtmp,r,p,t)
    Λ = gtmp
    @.. σ = sqrt(get_temperature(p, t)) * sqrtmass

    @views for i in axes(r, 3)
        @.. noise = σ * W.dW[:,:,i][:] / sqdt

        γ, c = LAPACK.syev!('V', 'U', Λ[:,:,i]) # symmetric eigen
        clamp!(γ, 0, Inf)
        for (j,i) in enumerate(diagind(c1))
            c1[i] = exp(-γ[j]*dt)
            c2[i] = sqrt(1 - c1[i]^2)
        end

        copyto!(flatvtmp, v[:,:,i])
        mul!(tmp1, transpose(c), flatvtmp)
        mul!(tmp2, c1, tmp1)
        mul!(flatvtmp, c, tmp2)
        
        mul!(tmp1, transpose(c), noise)
        mul!(tmp2, c2, tmp1)
        mul!(tmp1, c, tmp2)

        @.. flatvtmp += tmp1
        copyto!(v[:,:,i], flatvtmp)
    end

    RingPolymerArrays.transform_to_normal_modes!(r, p.beads.transformation)
    RingPolymerArrays.transform_to_normal_modes!(v, p.beads.transformation)
end

function step_O!(friction::LangevinCache, integrator, v::RingPolymerArray, r::RingPolymerArray, t)
    @unpack W, p, dt, sqdt = integrator
    @unpack c1, c2, sqrtmass, σ = friction

    @.. σ = sqrt(get_ring_polymer_temperature(p, t)) * sqrtmass

    for i in axes(r, 3)
        for j in RingPolymerArrays.quantumindices(v)
            @. v[:,j,i] = c1[i] * v[:,j,i] + c2[i] * σ[:,j] * W.dW[:,j,i] / sqdt
        end
    end

    for j in RingPolymerArrays.classicalindices(v)
        @. v[:,j,1] = c1[1] * v[:,j,1] + c2[1] * σ[:,j] * W.dW[:,j,1] / sqdt
    end
end


