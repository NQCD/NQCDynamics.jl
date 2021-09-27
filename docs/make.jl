using Documenter
using NonadiabaticMolecularDynamics
using NonadiabaticModels
using CubeLDFAModel

DocMeta.setdocmeta!(NonadiabaticMolecularDynamics, :DocTestSetup, :(using NonadiabaticMolecularDynamics); recursive=true)
DocMeta.setdocmeta!(NonadiabaticModels, :DocTestSetup, :(using NonadiabaticModels, Symbolics); recursive=true)

# Fix plots bug https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988
ENV["GKSwstype"] = "100"

@time makedocs(
    modules = [NonadiabaticMolecularDynamics, NonadiabaticModels, CubeLDFAModel],
    format = Documenter.HTML(
        prettyurls = true,
        canonical = "https://nqcd.github.io/NonadiabaticMolecularDynamics.jl/stable/",
        assets = ["assets/custom.css", "assets/favicon.ico"],
        ansicolor = true,
        ),
    clean = false,
    sitename = "NonadiabaticMolecularDynamics.jl",
    authors = "James Gardner and contributors.",
    pages = [
        "index.md"
        "getting_started.md"
        "NonadiabaticModels.jl" => Any[
            "nonadiabaticmodels/overview.md"
            "nonadiabaticmodels/analyticmodels.md"
            "Extra models and interfaces" => map(
                s -> "nonadiabaticmodels/models/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/nonadiabaticmodels/models")))
            )
        ]
        "Initial conditions" => Any[
            "initialconditions/dynamicaldistribution.md"
            "Sampling methods" => map(
                s -> "initialconditions/samplingmethods/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/initialconditions/samplingmethods")))
            )
        ]
        "Dynamics simulations" => Any[
            "dynamicssimulations/dynamicssimulations.md"
            "Dynamics methods" => map(
                s -> "dynamicssimulations/dynamicsmethods/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/dynamicssimulations/dynamicsmethods")))
            )
        ]
        "ensemble_simulations.md"
        "Examples" => [
            "examples/mdef.md"
            "examples/fssh.md"
            "examples/nrpmd.md"
            "examples/quantisation.md"
        ]
        "Developer documentation" => [
            "devdocs/new_methods.md"
            "devdocs/models.md"
        ]
        "API" => Any[
            "NonadiabaticModels" => map(
                s -> "api/nonadiabaticmodels/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/api/nonadiabaticmodels")))
            ),
            "NonadiabaticMolecularDynamics" => map(
                s -> "api/nonadiabaticmoleculardynamics/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/api/nonadiabaticmoleculardynamics")))
            ),
        ]
    ])

deploydocs(
    repo = "github.com/NQCD/NonadiabaticMolecularDynamics.jl",
)
