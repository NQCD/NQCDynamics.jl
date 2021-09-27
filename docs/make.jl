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
        "Dynamics simulations" => Any[
            "dynamicssimulations/dynamicssimulations.md",
            "Dynamics methods" => map(
                s -> "dynamicssimulations/dynamicsmethods/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/dynamicssimulations/dynamicsmethods")))
            )
        ]
        "initial_conditions.md"
        "ensemble_simulations.md"
        "Models" => [
            "models/overview.md"
            "models/model_library.md"
            "models/ldfa.md"
        ]
        "Examples" => [
            "examples/mdef.md"
            "examples/fssh.md"
            "examples/nrpmd.md"
            "examples/quantisation.md"
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
