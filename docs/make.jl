using Documenter, NonadiabaticMolecularDynamics
import JuLIP
using PyCall

doctestsetup = quote
    using NonadiabaticMolecularDynamics
end 
DocMeta.setdocmeta!(NonadiabaticMolecularDynamics, :DocTestSetup, doctestsetup; recursive=true)

@time makedocs(sitename="NonadiabaticMolecularDynamics.jl",
    modules=[NonadiabaticMolecularDynamics],
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
