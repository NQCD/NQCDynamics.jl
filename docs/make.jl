using Documenter, NonadiabaticMolecularDynamics

makedocs(sitename="NonadiabaticMolecularDynamics.jl",
    format=Documenter.HTML(prettyurls=true),
    pages=[
        "index.md"
        "getting_started.md"
        "initial_conditions.md"
        "ensemble_simulations.md"
        "calculators.md"
        "Models" => [
            "models/overview.md"
        ]
        "Dynamics" => [
            "dynamics/overview.md"
            "dynamics/classical.md"
            "dynamics/mdef.md"
        ]
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
