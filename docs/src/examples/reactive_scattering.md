# [Reactive scattering from a metal surface](@id example-h2scattering)

The current version of our Julia implementation allow us to simulate the state-to-state
vibrational de-excitation probability on reactive scattering events at metal surfaces for any diatomic molecule 
with a proper model to describe energies and forces (and eventually friction coefficients for MDEF simulations). 
Here, specifically, we examine the hydrogen reactive scattering on Ag(111) metal surface as a prototypical example to show how the general
workflow works.

To run this kind of simulations, a set of initial positions and velocities (``r_{0}^{i}`` and ``v_{0}^{i}``) with a specific
set of ro-vibrational quantum states ``\nu`` and ``j`` have to be initially generated
(see initial condition section/ sampling methods/ semiclassical EBK quantisation).
With a specific ro-vibrational quantum state is possible to compute different properties
after molecular collision and energy transfer with the metal surface like the vibrational
de-excitation probabilities discussed here.

In order to reproduce the state-to-state vibrational de-excitation probability results presented originally by [Maurer2019](@cite) for this system, the same initial conditions were generated through [`QuantisedDiatomic.generate_configurations`](@ref InitialConditions.QuantisedDiatomic.generate_configurations)
function with initial ro-vibrational quantum state (``\nu``=2 and ``j``=0) as was explored in the original paper. 

As shown earlier in the [EBK documentation](@ref ebk-sampling) we are able to generate
a semiclassically quantised distribution for a diatomic molecule on a collision course
with a metal surface.
Here we can follow that example with the [`H2AgModel`](@ref NNInterfaces.H2AgModel) model
to prepare our initial distribution and run our simulation.

Specifically, we have produced a set of initial conditions with different translational energy (`translational_energy` keyword) ranging from 0.2 to 1.4 eV, locating the hydrogen molecule at 8 Å away from the metal
surface (`height` keyword) with a normal incidence. Note, we use the unit transformation from Angstrom to atomic unit inside the function(`ang_to_au()`)

```@example h2scatter
using NQCDynamics
using NNInterfaces
using Unitful
using NQCDynamics.InitialConditions: QuantisedDiatomic
using JLD2

atoms = Atoms([:H, :H])
model = H2AgModel()
cell = PeriodicCell([11.1175 -5.5588 0.0; 0.0 9.628 0.0; 0.0 0.0 70.3079])
sim = Simulation(atoms, model; cell=cell)

ν, J = 2, 0           # selected ro-vibrational quantum states  
nsamples = 300        # number of configurations      
Ek = 0.5              # Translational energy [eV] ; range considered [0.2-1.4] eV
z = 8.0               # Height [Å]  ; fixed at 8 Å

configurations = QuantisedDiatomic.generate_configurations(sim, ν, J;
    samples=nsamples, translational_energy=Ek*u"eV", height= ang_to_au(z))
v = first.(configurations)
r = last.(configurations)

distribution = DynamicalDistribution(v, r, (3,2))
JLD2.save("distribution.jld2",Dict("dist"=>distribution,"atoms"=>atoms,"cell"=>cell))

```
The as-produced initial conditions can be subsequently saved in an external binary file
`distribution.jld2` format which store the distribution data (``r_{0}^{i}`` and ``v_{0}^{i}``) along
with useful information of the `unit cell` and `atoms type` with a dictionary structure. The
”jld2” format is a binary format which is able to store complex Julia data structure and it
is being widely used.

In order to produce unweighted distribution, the initial lateral and angular orientation were randomly selected within the unit cell. As example of the spacial and orientation distribution generated with this module,  a small subset of data (300 points) is shown in next Figure. To run our production simulations, however, a set of 80,000 initial velocities and positions were generated to save it to run ensemble simulation with the same pair of quantum numbers.

![initial conditions](../assets/figures/icond_scatter.png)

## MDEF with the LDFA

Now that we have set up the initial distribution and some of our simulation parameters,
we can choose which form of friction we would like use.
First, let's use the LDFA implementation provided by the
[CubeLDFAModel](@ref models-cubeldfa).
This takes a `.cube` file containing the electron density and will provide the friction
during the dynamics.
Here we have given the new model our model from above, which will provide the forces.

```@example h2scatter
using CubeLDFAModel
model = LDFAModel(model, "../assets/friction/test.cube", atoms, friction_atoms=[1,2], cell)
```

Now we can pass all the variables defined so far to the `Simulation` and run multiple
trajectories using [`run_ensemble`](@ref).

```@example h2scatter
sim = Simulation{MDEF}(atoms, model, cell=cell, temperature=300u"K")
ensemble = run_ensemble(sim, tspan, distribution;selection=1:20,
    dt=0.1u"fs", output=:position, trajectories=20, callback=terminate)
```

A full description of our reactive system also includes other basic variables definition
like type of atoms, unit cell and electronic temperature which can be properly defined
before to run any simulation. The generated distribution can be uploaded and use it in the following way.

```@example reading
system = load("distribution.jld2")

Tel = 0u"K"
cell = system["cell"]
distribution = system["dist"]
atoms = system["atoms"]
```

## MDEF with neural network friction 

Above, we used the LDFA interpretation of MDEF to perform the simulation.
However, the [`H2AgModel`](@ref NNInterfaces.H2AgModel) actually provides it's own
friction tensor trained on *ab initio* data.
This can be used by simply using the model directly, without wrapping it with the
[`LDFAModel`](@ref CubeLDFAModel.LDFAModel).

```@example h2scatter
model = H2AgModel()
sim = Simulation{MDEF}(atoms, model, cell=cell, temperature=300u"K")
ensemble = run_ensemble(sim, tspan, distribution;selection=1:20,
    dt=0.1u"fs", output=:position, trajectories=20, callback=terminate)
```

## Data analysis and truncation function

Since we are interested in the dynamics only when the molecule is close to the surface,
we can use a callback to terminate the simulation early to save us some time.
This requires defining a function that returns `true` when we want the simulation to
terminate.
This means we can set our time span relatively long since we expect most simulations to
terminate before reaching the time limit.

```@example h2scatter
using Statistics: mean
using LinearAlgebra: norm

h2distance(p) = norm(p[:,1] .- p[:,2])

function termination_condition(u, t, integrator)::Bool
    R = get_positions(u)
    zcom = au_to_ang(mean(R[3,:]))
    if zcom > 8.1
        return true
    elseif au_to_ang(h2distance(R)) > 2.5
        return true
    else
        return false
    end
end

terminate = DynamicsUtils.TerminatingCallback(termination_condition)
tspan = (0.0, 420.0u"fs")
nothing # hide
```

Specifically, for this example, the original criteria to define a reactive or scattering events was used and a
”termination condition function” was implemented as above. We consider a reactive event if the H-H
bond length is larger than 2.5 Å in any point of the molecular dynamics trajectory and a
scattering event if the ``H_2`` molecule scattering back with a vertical distance from the metal
surface larger than 8.1 Å.

As example to show this, we have run 20 trajectories with and without the truncation
function starting with  an initial translation energy at 1.0 eV. For both figures, the total and kinetic energies are shown in
the top panels along with the H-H distance and z coordinate of center of mass for each
individual trajectory.

![truncation](../assets/figures/scattering_truncation.png)



