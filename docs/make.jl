using Documenter
using DocumenterCitations
using NonadiabaticMolecularDynamics
using NonadiabaticModels
using CubeLDFAModel

DocMeta.setdocmeta!(NonadiabaticMolecularDynamics, :DocTestSetup, :(using NonadiabaticMolecularDynamics); recursive=true)
DocMeta.setdocmeta!(NonadiabaticModels, :DocTestSetup, :(using NonadiabaticModels, Symbolics); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"), sorting=:nyt)

@time makedocs(
    bib,
    sitename = "NonadiabaticMolecularDynamics.jl",
    modules = [NonadiabaticMolecularDynamics, NonadiabaticModels, CubeLDFAModel],
    strict = false,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://nqcd.github.io/NonadiabaticMolecularDynamics.jl/stable/",
        assets = ["src/assets/custom.css", "src/assets/favicon.ico"],
        ansicolor = true,
        ),
    authors = "James Gardner and contributors.",
    pages = [
        "Introduction" => "index.md"
        "Getting started" => "getting_started.md"
        "Atoms" => "atoms.md"
        "NonadiabaticModels.jl" => Any[
            "nonadiabaticmodels/overview.md"
            "nonadiabaticmodels/analyticmodels.md"
            "nonadiabaticmodels/ase.md"
            "nonadiabaticmodels/neuralnetworkmodels.md"
            "nonadiabaticmodels/frictionmodels.md"
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
        "Ensemble simulations" => "ensemble_simulations.md"
        "Examples" => map(
            s -> "examples/$(s)",
            sort(readdir(joinpath(@__DIR__, "src/examples")))
        )
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
        "References" => "references.md"
    ])


if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/NQCD/NonadiabaticMolecularDynamics.jl",
        push_preview=true
    )
end
