
export ang_to_au
export au_to_ang

export eV_to_au
export au_to_eV

export eV_per_ang_to_au
export au_to_eV_per_ang

export u_to_au
export au_to_u

export ps_inv_to_au
export au_to_ps_inv

export ps_to_au
export au_to_ps

export austrip_kwargs

auconvertstrip(u::Unitful.Units, x::Number) = ustrip(auconvert(u, x))

ang_to_au(x) = austrip(x*u"Å")
au_to_ang(x) = auconvertstrip(u"Å", x)

eV_to_au(x) = austrip(x*u"eV")
au_to_eV(x) = auconvertstrip(u"eV", x)

eV_per_ang_to_au(x) = austrip(x*u"eV/Å")
au_to_eV_per_ang(x) = auconvertstrip(u"eV/Å", x)

u_to_au(x) = austrip(x*u"u")
au_to_u(x) = auconvertstrip(u"u", x)

ps_inv_to_au(x) = austrip(x*u"ps^-1")
au_to_ps_inv(x) = auconvertstrip(u"ps^-1", x)

ps_to_au(x) = austrip(x*u"ps")
au_to_ps(x) = auconvertstrip(u"ps", x)

function austrip_kwargs(;kwargs...)
    stripped_kwargs = Dict{Symbol,Any}()
    for (key, val) in kwargs
        if val isa Unitful.Quantity
            stripped_kwargs[key] = austrip(val)
        else
            stripped_kwargs[key] = val
        end
    end
    stripped_kwargs
end
