using Documenter
using DocumenterCitations
using NonadiabaticMolecularDynamics
using NonadiabaticModels
using CubeLDFAModel, NNInterfaces

DocMeta.setdocmeta!(NonadiabaticMolecularDynamics, :DocTestSetup, :(using NonadiabaticMolecularDynamics); recursive=true)
DocMeta.setdocmeta!(NonadiabaticModels, :DocTestSetup, :(using NonadiabaticModels, Symbolics); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"), sorting=:nyt)

function find_all_files(directory)
    map(
        s -> joinpath(directory, s),
        sort(readdir(joinpath(@__DIR__, "src", directory)))
    )
end

@time makedocs(
    bib,
    sitename = "NonadiabaticMolecularDynamics.jl",
    modules = [NonadiabaticMolecularDynamics, NonadiabaticModels, CubeLDFAModel],
    strict = false,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://nqcd.github.io/NonadiabaticMolecularDynamics.jl/stable/",
        assets = ["assets/favicon.ico"],
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
            find_all_files("initialconditions/samplingmethods")
        ]
        "Dynamics simulations" => Any[
            "dynamicssimulations/dynamicssimulations.md"
            find_all_files("dynamicssimulations/dynamicsmethods")
        ]
        "Ensemble simulations" => "ensemble_simulations.md"
        "Examples" => find_all_files("examples")
        "Developer documentation" => find_all_files("devdocs")
        "API" => Any[
            "NonadiabaticModels" => find_all_files("api/nonadiabaticmodels")
            "NonadiabaticMolecularDynamics" => find_all_files("api/nonadiabaticmoleculardynamics")
        ]
        "References" => "references.md"
    ])


if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/NQCD/NonadiabaticMolecularDynamics.jl",
        push_preview=true
    )
end
