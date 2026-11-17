using Distributions
using Interpolations
using QuadGK
using NQCModels
using CairoMakie
using Unitful, UnitfulAtomic
using HokseonPlots, ColorSchemes, Colors, Printf
colorscheme = ColorScheme(parse.(Colorant, ["#045275", "#089099", "#7CCBA2", "#FCDE9C", "#F0746E", "#DC3977", "#7C1D6F"]))
colormap = HokseonPlots.NICECOLORS

fig = Figure(figure_padding=(1, 2, 1, 4), fonts=(;regular=projectdir("fonts", "MinionPro-Capt.otf")), size=(HokseonPlots.RESOLUTION[1]*3, HokseonPlots.RESOLUTION[2]*3))
ax1 = MyAxis(fig[1,1]; limits=(-27, 27, -1, 30), xgridvisible=false,ygridvisible=false)
ax2 = MyAxis(fig[2,1]; xlabel="Energy / eV", limits=(-27, 27, -1, 30), xgridvisible=false,ygridvisible=false)
axes = [ax1, ax2]
hidexdecorations!(ax1, ticks=false, minorticks=false)
sd = 0.2

function calculate_band_extrema(bath_dict::Dict{String,Any})
    @unpack centre, width = bath_dict
    bandmin = austrip(-(width/2 - centre) * u"eV")
    bandmax = austrip((width/2 + centre) * u"eV")
    return bandmin, bandmax
end

bath_dicts = Dict{String,Any}(
    "bandgap" => 5,
    "centre" => 0.0,
    "discretisation" => [:GapGaussLegendre, :GapTrapezoidalRule],
    "nstates" => 200,
    "width" => 50,
)

dicts = dict_list(bath_dicts)

for (i, bath_dict) in enumerate(dicts)

    @unpack bandgap, centre, discretisation, nstates, width = bath_dict

    bandmin = -(width/2 - centre) # unit is eV

    bandmax = (width/2 + centre) # unit is eV

    bandmin_val, bandmax_val = calculate_band_extrema(bath_dict)

    bath = eval(discretisation)(nstates, bandmin_val, bandmax_val, austrip(bandgap * u"eV"))

    Distribution_y = [Normal(center, sd) for center in ustrip(auconvert.(u"eV", bath.bathstates))]

    # Define a range of x values to evaluate the distributions
    DOS_x = collect(bandmin-20:0.01:bandmax+20)


    Coupling_y = ustrip(auconvert.(u"eV", bath.bathcoupling)).^2

    if length(bath.bathcoupling)==1
        Coupling_y = Coupling_y .* ones(length(Distribution_y))
    end 

    # Evaluate each Gaussian distribution at each x value and sum the results
    DOS_y = sum([pdf.(Distribution_y[i], DOS_x) .* Coupling_y[i] for i in eachindex(Distribution_y)])
    
    lines!(axes[i], DOS_x, DOS_y , color=colorscheme[i], label= string(discretisation))
    poly!(axes[i], DOS_x, DOS_y , color=(colorscheme[i], 0.3), transparency=true)

    axislegend(axes[i], position=:rt)
    
end
Label(fig[:,0], "Density / arb. u", rotation=π/2, padding=(0, 0, 0, 0))

save(projectdir("plots/DOS_bath", "DOS_bath_discretisation_compare_test.svg"),fig)

fig