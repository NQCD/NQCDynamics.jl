```@setup logging
@info "Expanding src/dynamicssimulations/dynamicsmethods/iesh.md..."
start_time = time()
```
# [Independent electron surface hopping (IESH)](@id iesh-dynamics)

Independent electron surface hopping (IESH) ([Shenvi2009](@cite)[Gardner2023](@cite)) is an extension of Tully's fewest-switches surface hopping (FSSH) [Tully1990](@cite) that describes a single particle (single electronic state) interacting with a metal, represented as a bath of electronic states.

IESH is classified as a mixed-quantum-classical method, where the nuclei are treated classically and the electrons are treated quantum mechanically. Like FSSH and other such surface hopping methods, the nuclei evolve on a single adiabatic potential energy surface at any given moment. At each timestep, a hopping probability is evaluated. If the hopping probability is larger than a uniform random number between 0 and 1, the active state is switched and the adiabatic propagation continues on the new electronic state.
When this algorithm is applied to an ensemble of trajectories, the discrete adiabatic state populations approximate the quantum mechanical populations for each state.

The key difference with IESH comes from the assumption that electrons in the systems are independent. This allows:
 1. the many-electron Hamiltonian of the system to be decomposed into a sum of single electron Hamiltonians ``\hat{H}_{el}^{1}(\mathbf{q})``,
 2. the many-electron wavefunction ``\ket{\Psi}`` can be fomrulated as a single slater determinant ``\ket{\psi_{1} \psi_{2} ... \psi_{N_{e}}}``


The IESH classical Hamiltonian can be written as
```math
H_{\textrm{IESH}}(\mathbf{q}(t), \mathbf{p}) = \sum_{\alpha=1}^{N} \frac{p_{\alpha}^{2}}{2m_{\alpha}} + U_{0}(\mathbf{q}_{\alpha}(t)) + \sum_{k \in \mathbf{s}(t)} \lambda_{k} (\mathbf{q}_{\alpha}(t))
```
where ``\mathbf{p}`` is the vector of momenta, ``\mathbf{q}`` the positions and ``|\alpha`` denotes the nuclear degree of freedom.
``\mathbf{s}(t)`` is a binary vector that indicates the currently occupied adiabatic electronic states as a function of time
As such, this Hamiltonian describes classical dynamics that proceeds under the influence
of the model's ground state potential ``U_{0}(\mathbf{q}_{\alpha}(t))``, adjusted by the sum of the eigenvalues of ``\hat{H}_{el}^{1}(\mathbf{q})`` to corresponding to occupied electronic states given in ``s(t)``.

Similarly to FSSH, ``s(t)`` must be obtained as the dynamics progress since it itself is time dependent. Again, this is achieved stochastically for each trajectory by making probabilistic hops between surfaces.
The probabilities for making successful hops are obtained by integrating the electronic Schrödinger equation alongside the dynamics as
```math
i \hbar \dot{c}_{k}^{l}(t) = \lambda_{k}(\mathbf{q}(t))c_{k}^{l}(t) - i \hbar \sum_{\alpha=1}^{N}\sum_{j \neq k} \frac{p_{\alpha}}{m_{\alpha}} d_{\alpha,jk}(\mathbf{q}(t)) c_{j}^{l}(t)
```
In this equation, ``c_{k}^{l}(t)`` are the time dependent expansion coefficients for state ``k`` and electron ``l``.
``\mathbf{d}_{\alpha,jk}`` is the nonadiabatic coupling between adiabatic states ``j`` and ``k``.
The hopping probability (which corresponds to replacing index ``k`` in ``s(t)`` with ``j``) is calculated as
```math
g_{k \to j} = \textrm{max}\left( \frac{B_{jk} \Delta t}{A_{kk}}, 0 \right)
```
Where
```math
\begin{align*}
    B_{jk} &= -2 \textrm{Re} (A_{kj}^{*}) \sum_{\alpha = 1}^{N} \frac{p_{\alpha}}{m_{\alpha}} d_{\alpha,jk}\\
    A_{kj} &= \braket{\mathbf{k}|\psi} \braket{\psi | \mathbf{j}}  
\end{align*}
```
At each timestep, a uniform random number between 0 and 1 is generated which is compared to the calculated probabilities. If the probability is higher than the random number, then a hop is attempted.

In addition, after a hopping event is selected to occur, the velocity needs to be rescaled in the direction of the nonadiabatic coupling vector. A selected hop will only be successful when there is sufficient kinetic energy for the rescaling such that energy is conserved after the hop.

If there is insufficient kinetic energy, this is termed a frustrated hop, and the dynamics proceeds without performing a hop. When a hop is successful, the kinetic energy is adjusted and the occupation vector ``s(t + \Delta t)`` is modified to reflect the change in state occupation.

## Algorithm

1. Integrate classical dynamics for one timestep
2. Integrate electronic dynamics for one timestep
3. Evaluate hopping probability
4. Perform hop if sufficient probability and kinetic energy
5. Rescale velocity if hop is performed
6. Return to step 1

!!! note

    With `DifferentialEquations.jl` we use a callback to perform the surface hopping
    procedure such that steps 1 and 2 are performed by the DifferentialEquations solvers and steps 3, 4, 5 are
    performed by the callback.

More details of this method's implementation within NQCD, and associated simulation results for molecular scattering and thermal desorption from a metal surface using IESH were published at the following reference ([Gardner2023](@cite)).

## Interfacing with the `AndersonHolstein` model

The Newns-Anderson or Anderson-Holstein model ([Anderson1961](@cite)[Newns1969](@cite)) is a commonly used impurity model which has many applications to surface systems where one wants to simulate an absorbate molecule interacting with a surface. This model works well with the IESH method due to presence of an the explicit bath. 

As such, the IESH single electron Hamiltonian within the framework of the Newns-Anderson model can be written as
```math
    \hat{H}_{el}^{1}(\mathbf{q}) = h(\mathbf{q})\ket{d}\bra{d} + \sum_{k=1}^{N} \epsilon_{k} \ket{k}\bra{k} + \sum_{k=1}^{N} V_{k}(\mathbf{q})\left( \ket{k}\bra{d} + \ket{d}\bra{k} \right)
```
Where ``\ket{d}`` and ``\ket{k}`` are the single electron states, ``\ket{\psi_{l}}`` acted upon by the creation operators for an electron in the absorbate state (``\hat{d}``) and ``k^\textrm{th}`` bath state (``\hat{c}_{k}^{\dagger}``) respectively.
``h(\mathbf{q})`` is the difference in energy between the absorbate ground and excited states (``U_{1}(\mathbf{q}) - U_{0}\mathbf{q}``), ``\epsilon_{k}`` is the energy contribution from the ``k^\textrm{th}`` bath state, and ``V_{k}(\mathbf{q})`` is the coupling strength between the ``k^\textrm{th}`` bath state and the absorbate state.

More details of the Newns-Anderson model and it's implementation within NQCD is provided [here](../../NQCModels/systembathmodels.md).

## Example

In this section we can investigate the results for molecular scattering from a metal surface, obtained for a single trajectory using IESH.

First, the simulation parameters are created. Here, we have a single atom with a mass of
`10.54` a.u. and we are using `ErpenbeckThoss` model ([Erpenbeck2018](@cite)[Erpenbeck2019](@cite)), provided by [NQCModels.jl](@ref), to describe the absorbate, and the `AndersonHolstein` model to provide the coupling to a metallic surface. The Gauss-Legendre bath discretisation method (`ShenviGaussLegendre()`) implemented by Shenvi ([Shenvi2009](@cite)) has been used in preference to a uniform discretisation to better characterise the electronic states close to the fermi level. More details regarding discretisation for system-bath models can be found [here](../../NQCModels/systembathmodels.md).
```@example iesh
using Random; Random.seed!(10)
using NQCDynamics
using Unitful, UnitfulAtomic

atoms = Atoms(2000)

#### Absorbate
Γ = austrip(0.2u"eV") # Defines coupling strength between diabatic ground state and diabatic excited state of the absorbate
thossmodel = ErpenbeckThoss(;Γ, m=atoms.masses[1]) # This model defines the absorbate

#### Bath
nstates = 20 # Number of states into which the electronic continuum of the substrate is discretised into
bandwidth = austrip(50*u"eV") #Spectral width of the electronic bands
bandmin = -bandwidth / 2
bandmax = bandwidth / 2
bath = ShenviGaussLegendre(nstates, bandmin, bandmax) #Discretisation method

#### Combined Model
model = AndersonHolstein(thossmodel, bath)
```

As this is mixed-quantum-classical dynamics method where both the electrons and nuceli are evolved, initial conditions must be generated both for the nuclei (`nuclear_distribution`) and electrons (`electronic_distrbution`) and passed as a "product distribution" to the `run_dynamics()` call. 
With an absorbate-bath system, the electronic distribution is given by some Fermi-Dirac distribution generated at a given temperature, from which the electronic state occupations can be sampled during initialisation.
```@example iesh
####  nuclear parameters
r0 = austrip(5u"Å") # Projectile is 5A away from surface
m = atoms.masses[1] # Mass of projectile
ke = austrip(2u"eV") # Initial kinetic energy
v0 = -sqrt(2*ke/m) # Velocity direction towards the surface with magnitude defined by kinetic energy
nuclear_distribution = DynamicalDistribution(v0, r0, (1,1))

#### electronic parameters
fermi_level = 0*u"eV" # We reference the spectrum of electronic states of the surface to the Fermi level
temperature = 300u"K"
electronic_distribution = FermiDiracState(fermi_level, temperature)

#### product distribution
dist = electronic_distribution * nuclear_distribution
```

The final component is defining the "termination condition" that bring the simulation to an end. For this scattering simulation the run is terminated when the absorbate is more the 5.5 Å away from the surfce with a non-zero velocity, or, if the run time exceeds some cut-off time that is sensible.
```@example iesh
### Defining termination criterion
dcut = 3*r0
tcut = abs(dcut/v0)

function termination_condition(u, t, integrator)::Bool

    return ((t > tcut) || ((mean(DynamicsUtils.get_positions(u)) > austrip(5.5u"Å")) && (mean(DynamicsUtils.get_velocities(u)) > 0)))
end
terminate = DynamicsUtils.TerminatingCallback(termination_condition)
```

Finally, the trajectory can be run by passing all the parameters we have set up so far.
Here, we request both the `OutputDiscreteState` output which is equal to ``s(t)``, `OutputDiabaticPopulation`, which gives us the population of each diabatic state along the trajectory and `OutputSurfaceHops`, which details the number of hops between iesh surfaces were allowed during the trajectory.
```@example iesh
using Statistics
ntrajs = 1
dt = 0.01u"fs"
output= (OutputPosition, OutputVelocity, OutputTotalDiabaticPopulation, OutputDiabaticPopulation, OutputAdiabaticPopulation, OutputTotalAdiabaticPopulation, OutputSurfaceHops)

### Run IESH simulations
sim = Simulation{AdiabaticIESH}(atoms, model)

kick = run_dynamics(sim,
                    (0.0, dt),
                    dist;
                    trajectories=1,
                    callback=terminate,
                    output=output,
                    dt=dt,
                   )

result = run_dynamics(sim,
                    (0.0, tcut),
                    dist;
                    trajectories=ntrajs,
                    callback=terminate,
                    output=output,
                    saveat=[0,tcut],
                    dt=dt,
                   )
result
```

Now we can plot ``s(t)`` throughout the trajectory, and read out the number of hops that were allowed.
```@example iesh
println("s(t) = $(result[:OutputDiscreteState])")
println("Num. surface hops = $(result[:OutputSurfaceHops])")
```

Similarly, we can plot the diabatic populations. Since FSSH is performed in the adiabatic
representation, even in the case of few hops, the diabatic populations can look dramatically
different depending on the complexity of the model Hamiltonian. 
```@example iesh
using Plots
plot(result[:OutputDiabaticPopulation][1], xlabel="bath states / eV", ylabel="Diabatic Population", label="start")
plot!(result[:OutputDiabaticPopulation][end], label="end")
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
