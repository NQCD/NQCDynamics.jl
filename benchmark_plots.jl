import JSON
using CairoMakie
using Glob

# Get information from outside the plotting script: What is the newest version and what are we comparing to?
newest_version = get(ENV, "CURRENT_VERSION", "v0.15.1") # A known tested version
github_context = get(ENV, "NEW_VERSION", "new") # The new thing to test. 
graphics_path = get(ENV, "GRAPHICS_PATH", "benchmark-graphics") # The new thing to test. 

# Base path of the benchmark repository
benchmark_path = ENV["BENCHMARK_REPO"]

# Find all JSON files for current and new versions
newest_path = "$(benchmark_path)/NQCDynamics.jl/commits/$(github_context)/"
current_path = "$(benchmark_path)/NQCDynamics.jl/tags/$(newest_version)/"

if isdir(newest_path)
	newest_files = split.(Glob.glob("*.json", "$(newest_path)"), "/") .|> last
else
	newest_files = String[]
end
if isdir(current_path)
	current_files = split.(Glob.glob("*.json", "$(current_path)"), "/") .|> last
else
	current_files = String[]
end

println(newest_files)

base_theme = Theme(
  palette = (
	marker = [:circle, :rect, :utriangle, :cross, :star5],
	linestyle = [:solid, (:dot, :dense), (:dash, :dense), (:dashdot, :dense), (:dashdotdot, :dense)],
  ),
  CairoMakie=(px_per_unit = 10.0, ), # Store raster images with sufficient resolution
  size = 28.1 .* (12,9),
  Figure = (
	size = 28.1 .* (12,9), # 28.1 pt per cm, so this should be 12x8 cm
	fontsize=12, # Match PDF font size
	figure_padding = 0 # Let whatever is used to display images do the padding. 
  ),
  Axis = (
	xgridvisible=false, # remove gridlines
	ygridvisible=false,
	xminorticksvisible=true, # minor ticks are nice
	yminorticksvisible=true,
	# ticks face out by default
	spinewidth = 0.75 # Slightly thinner line thickness
  ),
  Colorbar = (
	spinewidth = 0.75,
  ),
  Legend = (
	padding = 2, # reduce overall padding around Legends
	framevisible = false, # remove frame
	patchsize = (15,15),
	fontsize = 8,
  ),
  Errorbars = (
	whiskerwidth = 2.5,
	strokecolor = "black",
  ),
  Rangebars =(
	strokewidth = 0.5,
	strokecolor = "black",
  ),
  Density =(
	strokewidth = 0.75,
	strokearound = true,
  ),
  Lines = (
	linewidth = 1.5,
	cycle = Cycle([:color, :linestyle], covary = true)
  ),
  Scatter = ( # Give markers a black border
	markersize = 8,
	cycle = Cycle([:color, :marker], covary = true),
	strokewidth = 0.5,
	strokecolor = "black",
  ),
  Hist = ( # Histograms should have thin black borders
	strokewidth = 0.5,
	strokecolor = "black",
	cycle = :color
  ),
  BarPlot = ( # Bar plots should have thin black borders
	strokewidth = 0.5,
	strokecolor = "black",
	cycle = :color
  ),
  ScatterLines = (
	linestyle = :dash,
	linewidth = 0.75,
	markersize = 10,
	cycle = Cycle([:color, :marker, :linestyle], covary = true),
	strokewidth = 0.5,
	strokecolor = "black",
  ),
  Heatmap = (
	cycle = Cycle([:colormap]),
  ),
)

CairoMakie.update_theme!(base_theme)

# Navigate through each new version file
for filename in newest_files
	newest_dict = JSON.parsefile(newest_path*filename)
	# Check if an old version equivalent exists
	current_dict = !isnothing(findfirst(x->x==filename, current_files)) ? JSON.parsefile(current_path*filename) : Dict{String, Any}()
	
	labels = symdiff(union(keys(current_dict), keys(newest_dict)), ["title_for_plotting"]) |> collect
	
	times_new = [haskey(newest_dict, k) ? newest_dict[k]["Time"] : 0.0 for k in labels]
	times_current = [haskey(current_dict, k) ? current_dict[k]["Time"] : 0.0 for k in labels]
	allocs_new = [haskey(newest_dict, k) ? newest_dict[k]["Allocs"] : 0.0 for k in labels]
	allocs_current = [haskey(current_dict, k) ? current_dict[k]["Allocs"] : 0.0 for k in labels]
	
	# Remove infinity for newly added tests
	relative_performance_new = times_new ./ times_current
	relative_performance_new[findall(relative_performance_new .== Inf)] .= 1.0
	relative_performance_old = times_current ./ times_current
	
	#all_values = vcat(relative_performance_new, relative_performance_old)
	all_values = vcat(times_new, times_current)
	category = repeat(collect(1:length(labels)), 2)
	dodge = repeat([1,2], inner = length(relative_performance_new))
	#  Conditionally format performance changes red and green. 
	colors = [colorant"#25afdc" for i in 1:2*length(relative_performance_new)]
	# ≥5% slowdown - red
	colors[findall((relative_performance_new .- relative_performance_old) .≥ 0.05)] .= colorant"#e14343"
	# Minor slowdown or improvement - green
	colors[findall((relative_performance_new .- relative_performance_old) .≤ 0.05)] .= colorant"#1ccf14"
	# New functionality - always show green
	colors[findall(isnan.(relative_performance_new .- relative_performance_old))] .= colorant"#1ccf14"
	
	# Plot a histogram of test timings where the old time is equivalent to 100%. 
	figure = Figure(
		size = 28.1 .* (
			# any(length.(labels) .≥ 12) ? maximum(length, labels) : 12,
			max(all_values...) ≥ 12.0 ? max(all_values...) : 12,
			length(labels) > 2 ? 4*length(labels) : 8,
		),
	)
	Label(figure[1,1], text = newest_dict["title_for_plotting"],)
	benchmark = Axis(
		figure[3,1],
		xlabel = "Time / s",
		ylabel = "Test set",
		xautolimitmargin = (0.0,0.1),
		yautolimitmargin = (0.15,0.15),
		yticks = (1:length(labels), join.(split.(labels), "\n")),
		# xticks = (0:0.25:1.0, string.(collect(0:25:100))),
		yminorticksvisible = false,
	)
	Legend(
		figure[2,1],
		[
			PolyElement(color = colorant"#25afdc", strokecolor = colorant"black", strokewidth=0.75), 
			[PolyElement(color = colorant"#e14343", strokecolor = colorant"black", strokewidth=0.75, points=Point2f[(0, 0), (1, 0), (0, 1)]), PolyElement(color = colorant"#1ccf14", strokecolor = colorant"black", strokewidth=0.75, points=Point2f[(1, 1), (1, 0), (0, 1)])]
		],
		[newest_version, length(github_context) > 7 ? github_context[1:7] : github_context],
		height = 25,
		orientation = :horizontal,
	)
	bp = barplot!(
		benchmark,
		category,
		all_values,
		dodge = dodge,
		color = colors,
		direction = :x,
		bar_labels = ["$(round(t; digits=3)) s" for t in Iterators.flatten([times_new, times_current])],
		strokecolor = colorant"black",
		strokewidth = 0.5,
		flip_labels_at = 0.75 * max(all_values...),
	)
	rowgap!(figure.layout, 2)
	colsize!(figure.layout, 1, Relative(1.0))
	
	save("$(graphics_path)/$(newest_dict["title_for_plotting"]).png", figure)
	
	# Output warnings for files present in old version that don't have a new one. 
end


