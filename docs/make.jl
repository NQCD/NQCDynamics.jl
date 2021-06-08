using Documenter, NonadiabaticMolecularDynamics
using NonadiabaticModels
import CubeLDFAModel
using PyCall
using JuLIP

doctestsetup = quote
    using NonadiabaticMolecularDynamics
    using NonadiabaticModels
end 
DocMeta.setdocmeta!(NonadiabaticMolecularDynamics, :DocTestSetup, doctestsetup; recursive=true)
DocMeta.setdocmeta!(NonadiabaticModels, :DocTestSetup, :(using NonadiabaticModels); recursive=true)

@time makedocs(sitename="NonadiabaticMolecularDynamics.jl",
    modules=[NonadiabaticMolecularDynamics, NonadiabaticModels, CubeLDFAModel],
    format=Documenter.HTML(prettyurls=true, assets=["assets/custom.css"]),
    pages=[
        "index.md"
        "getting_started.md"
        "initial_conditions.md"
        "ensemble_simulations.md"
        # "calculators.md"
        "Models" => [
            "models/overview.md"
            "models/model_library.md"
            "models/ldfa.md"
        ]
        # "Dynamics" => [
        #     "dynamics/overview.md"
        #     "dynamics/classical.md"
        #     "dynamics/mdef.md"
        # ]
        "Examples" => [
            "examples/mdef.md"
            "examples/fssh.md"
            "examples/nrpmd.md"
            "examples/quantisation.md"
        ]
    ])

deploydocs(
    repo = "github.com/maurergroup/NonadiabaticMolecularDynamics.jl",
)
