using Documenter, NonadiabaticMolecularDynamics

makedocs(sitename="NonadiabaticMolecularDynamics.jl",
    format=Documenter.HTML(prettyurls=false),
    pages=[
        "NonadiabaticMolecularDynamics.jl Documentation" => "index.md"
        "calculators.md"
        "Models" => [
            "models/overview.md"
        ]
        "Dynamics" => [
            "dynamics/overview.md"
            "dynamics/classical.md"
            "dynamics/mdef.md"
        ]
        "Initial Conditions" => [
            "initial_conditions/initial_conditions.md"
        ]
    ])

deploydocs(
    repo = "github.com/maurergroup/NonadiabaticMolecularDynamics.jl",
)
