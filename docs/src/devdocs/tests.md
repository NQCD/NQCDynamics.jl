## Unit test overview

This page lists all unit tests in `@test/`, grouped by category and file. Each bullet corresponds to an `@testset` name (nested testsets are shown indented). Use this overview to quickly locate relevant tests when developing or debugging.

- How to run locally:
  - `julia --project=test -e 'using Pkg; Pkg.test()'`
  - Or from the project: `julia --project -e 'using Pkg; Pkg.test()'`

### InitialConditions

- `test/InitialConditions/advancedmh_sampling.jl`
  - Classical
    - get_proposal — Checks proposal generation for classical MH/HMC samplers produces finite, well-shaped moves consistent with target distributions.
    - Energy expectation — Verifies sampled kinetic/potential energy expectations match analytical results for classical ensembles.
    - HMC gradient — Confirms force/gradient usage in HMC proposals is consistent and stable; catches sign or scaling errors.
    - Energy expectation HMC — Ensures HMC sampling recovers correct energy moments, indicating proper integrator and acceptance logic.
  - Ring polymer
    - get_proposal — Validates path-integral proposals (centroid/internal modes) preserve correct bead correlations for PIMD.
    - Energy expectation — Checks ring-polymer energy moments at finite T align with quantum-statistical expectations [@Craig2004; @Korol2020].
  - Free ring polymer radius_of_gyration — Confirms free-RP scaling of radius of gyration vs temperature/beads behaves as expected for a harmonic free chain [@Tuckerman2010].

- `test/InitialConditions/monte_carlo.jl`
  - propose_move! — Proposal mechanism produces symmetric, bounded perturbations with correct dimensionality.
  - assess_proposal! — Metropolis–Hastings acceptance decision computed from energy differences with proper temperature scaling.
  - write_output! — Output writer produces structured trajectories/observables required for downstream analysis.
  - acceptance_probability — Analytical acceptance probability limits tested on simple ΔE values.
  - run_monte_carlo_sampling — End-to-end MC sampling generates distributions with stable acceptance rates and equilibrated energies.
  - run_monte_carlo_sampling — Path-integral MC variant behaves analogously but over bead-expanded state.
  - propose_centroid_move! — Centroid-only proposals affect only centroid modes as required for constrained/path-integral updates.

- `test/InitialConditions/quantised_diatomic.jl`
  - separate/combine_slab_and_molecule — Asserts index bookkeeping correctly splits and recombines molecule and slab atoms in 1D/3D cases.
    - 3D, with slab — Honors slab indices and ordering.
    - 3D, no slab — Molecule-only geometry consistent.
    - 1D, no slab — Reduced DoF indices handled correctly.
  - calculate_force_constant — Extracts curvature from local potential; validates against harmonic limit.
  - subtract_centre_of_mass! — Removes overall COM translation while preserving internal coords.
  - apply_random_rotation! — Applies uniform random rotations; preserves bond lengths/energies.
  - generate and check results — Generates quantized vibrational configs and checks internal consistency [EBK approach: @Larkoski2006].
  - position_above_surface! — Places molecule above slab at specified height without overlap.
  - velocity_from_energy — Converts scalar kinetic energy to velocity magnitude respecting masses.
  - apply_translational_impulse! — Impulses adjust linear momentum without changing internal modes.
  - find_total_energy — Sums kinetic and potential energies; used by acceptance checks.
  - find_integral_bounds — Determines integral bounds for vibrational states from turning points.
  - calculate_diatomic_energy
    - 3D, no slab — Vibrational/rotational energy components combine to expected totals at given geometry.
  - assemble_evaluation_geometry
    - 3D, no slab — Assembles geometry/environment container with consistent dimensions.
    - 1D, no slab — 1D specialization uses scalar axes consistently.
  - build_molecule
    - 3D, no slab — Constructs `Atoms` with correct species, positions, velocities.
    - 1D, no slab — 1D positions represented consistently with modeling code.
  - fit_binding_curve and plot_binding_curve — Fits Morse/harmonic parameters to sampled points and checks residuals are small.
  - generate_configurations 1D
    - forward and backward: ν = $quantum_number — Forward/backward sampling for vibrational quantum number ν yields symmetric distributions.

### Ensembles

- `test/Ensembles/ensembles.jl`
  - run_dynamics — End-to-end ensemble execution returns outputs of expected shape/types.
  - run_dynamics — Variant path checks scheduling and result collation.
  - run_dynamics, reduction=$reduction — Mean/Sum reductions produce correct aggregated observables.

- `test/Ensembles/reductions.jl`
  - $reduction — Verifies that reduction operators are associative and produce expected results on toy data.
  - FileReduction — Streaming reductions to files preserve values and order.

- `test/Ensembles/selections.jl`
  - OrderedSelection — Deterministic selection order for reproducibility.
  - RandomSelection — Stochastic selection yields unbiased coverage with fixed seed reproducibility.

- `test/Ensembles/outputs.jl`
  - OutputFinal — Collects final-time observables; validates schema/units.
  - OutputDissociation — Flags dissociation events based on distance/energy thresholds.
  - OutputQuantisedDiatomic — Outputs state populations for quantized diatomic prep.
  - PopulationCorrelationFunction — Computes state population auto/cross-correlations; sanity-checked against trivial cases.
    - Ehrenfest — Mean-field electronic population definition [@Tully1990].
    - FSSH — Hopping-based state populations [@Tully1990; @Subotnik2016; @landry2013].
    - eCMM — Mapping-based electronic populations [@He2019; @gao2020].

### Analysis

- `test/Analysis/diatomic.jl`
  - Desorption Analysis — Post-processing extracts desorption probabilities/energetics; checks consistency with trajectory flags.

### Structure

- `test/Structure/structure.jl`
  - Minimum distance translation — MIC shifts minimize interatomic distances under PBC.
  - PBC-compatible functions — PBC-safe distance and vector ops return correct values on edge cases.

### Core

- `test/Core/FastDeterminant.jl`
  - det! $T — Fast determinant routine returns correct values and types across real/complex element types.

- `test/Core/estimators.jl`
  - @estimate — Macro expansion captures the intended observable with no allocations beyond baseline.
  - potential_energy — Potential estimator matches model evaluations; checks units and shapes.
  - kinetic_energy — Kinetic estimator uses masses and velocities correctly in classical limit.

- `test/Core/simulations.jl`
  - Thermostats — BAOAB/GLA-type thermostats preserve correct stationary distribution and target temperature [@Leimkuhler2012; @Tuckerman2010].
  - get_temperature(Simulation) — Temperature estimator yields target T in equilibrium.
  - get_ring_polymer_temperature — Centroid/internal-mode temperatures reported correctly.
  - nfunctions and size — Bookkeeping: dimensionalities and function counts consistent with model and atoms.
  - masses — Mass vectors/tensors match atoms and bead expansion.
  - mobileatoms and distribution generation — Mobility masks and distribution factories generate consistent initial conditions.

- `test/Core/calculators.jl`
  - General constructors — Multiple dispatch yields correct calculator/cache types for model/simulation variants.
  - ClassicalModel_Cache — Classical single-bead cache produces energies/forces with zero extraneous allocations.
  - RingPolymer_ClassicalModel_Cache — RP cache aggregates beadwise quantities and exposes centroid projections.
  - Abstract_QuantumModel_Cache
    - Potential evaluation — Diabatic/adiabatic potentials consistent with model definitions.
    - Dependent evaluation — Derived quantities (e.g., gradients, couplings) update properly.
    - Zero allocations — Hot loops allocation-free to ensure performance stability.
  - RingPolymer_QuantumModel_Cache
    - Potential evaluation — Per-bead quantum potentials consistent across beads.
    - Dependent evaluation — NACs and derived data consistent per bead.
    - Centroid dependent evaluation — Centroid-projected observables computed correctly.
    - Extra quantities — Additional diagnostics present and consistent.
    - Zero allocations — Allocation-free property maintained in tight loops.
  - Abstract_QuantumModel_Cache — Cross-checks abstract interface invariants.
  - RingPolymer_QuantumModel_Cache — Same invariants for RP specialization.
  - Friction_Cache — Friction tensors and noise scales are shaped/symmetric and positive semidefinite where required [@Gerrits2020; @Box2023].
  - RingPolymer_Friction_Cache — RP friction expansion preserves centroid/internal mode coupling structure.

### Dynamics

- `test/Dynamics/DynamicsUtils.jl`
  - classical_potential_energy : $(name(sim)) — Energy helpers return scalar energies consistent with model calculators.
  - classical_kinetic_energy : $(name(sim)) — Kinetic energy computed from velocities/masses correctly.
  - divide_by_mass! : $(name(sim)) — Mass scaling of vectors works in-place with correct broadcasting.
  - transform_density! — Density transformation routines are unitary/stable.
  - initialise_adiabatic_density_matrix
    - Adiabatic — Initializes pure-state density matrices on selected adiabatic surface.
    - Diabatic — Initializes diabatic density consistent with chosen basis.

- `test/Dynamics/bcbwithtsit5.jl`
  - Ehrenfest — Fixed-step vs Tsit5 comparisons for reference trajectories; checks mean-field equations [@Tully1990].
  - FSSH — Hopping dynamics integrated by Tsit5 matches constraints and event handling [@Tully1990; @Subotnik2016].

- `test/Dynamics/bcme.jl`
  - BCME — Validates Broadened CME limiter and detailed balance in simple models; basic sanity versus CME.

- `test/Dynamics/cme.jl`
  - Phonon relaxation Tᵢ=$(T)T — Temperature-dependent relaxation towards equilibrium rates is monotonic and physically plausible.

- `test/Dynamics/cmm.jl`
  - generate_random_points_on_nsphere — Uniform sampling on hypersphere surface validated by moments.
  - Population correlation — Electronic population correlation normalized at t=0 and decays smoothly.
    - Gamma = $γ — Damping parameter γ controls decay rate as expected.

- `test/Dynamics/diabatic_mdef.jl`
  - Friction comparison — Diabatic MDEF friction consistent with analytic/LDFA-inspired references in limits [@Gerrits2020; @Maurer2019].
  - Simulation{DiabaticMDEF} — Problem construction wires friction/noise consistent with fluctuation–dissipation.
  - RingPolymerSimulation{DiabaticMDEF} — RP analogue preserves centroid temperature and bead noise scaling.

- `test/Dynamics/electronic_dynamics.jl`
  - ElectronicODEProblem, norm preservation — Wavefunction norm preserved under unitary evolution.
  - DensityMatrixODEProblem, norm preservation — Trace of density matrix conserved.

- `test/Dynamics/ehrenfest.jl`
  - Ehrenfest — Mean-field nuclei with time-dependent electronic expectation; sanity checks vs model predictions [@Tully1990].
  - get_diabatic_population — Correct transformation from adiabatic coefficients to diabatic populations.
  - Algorithm comparison — Compares integrators for stability/accuracy on same problem.
  - FermiDirac Ehrenfest algorithm comparison — Checks finite-temperature Fermi-Dirac occupancy handling is consistent.
  - Spin boson population dynamics — Reproduces expected population oscillations/relaxation trends of the spin-boson model.
  - Ehrenfest RPMD — Combines mean-field with RP modes; verifies conservation/stability in small systems [@Richardson2013; @Richardson2017].

- `test/Dynamics/ehrenfest_na.jl`
  - Algorithm comparison: $method — Cross-method consistency check in nonadiabatic mean-field limit.

- `test/Dynamics/fssh.jl`
  - FSSH — Validates core components of Tully’s algorithm [@Tully1990; @Subotnik2016].
    - get_diabatic_population — Population estimator matches coefficient bookkeeping [@landry2013].
    - select_new_state — Hops sampled from cumulative probabilities; edge cases behaved deterministically.
    - rescale_velocity! — Momentum rescaling conserves total energy along coupling direction; frustrated hops handled.
    - calculate_potential_energy_change — ΔE computed from surface energies for hop tests.
    - execute_hop! — State, momentum, and bookkeeping updated atomically and consistently.
    - run_dynamics — End-to-end hop dynamics produces stable trajectories and populations.
  - RPSH — Ring Polymer Surface Hopping components [@Shushkov2012; @Shakib2017].
    - rescale_velocity! — RP-rescaling impacts centroid/internal modes as specified in RPSH.
    - calculate_potential_energy_change — ΔE computed per bead and combined appropriately.
    - execute_hop! — RP state updates honor bead coupling.
    - run_dynamics — RPSH trajectories remain stable under test parameters.

- `test/Dynamics/iesh.jl`
  - DynamicsVariables — IESH state vector and caches initialized properly [@Shenvi2009; @Gardner2023].
  - set_unoccupied_states! — Unoccupied manifold set consistently with electron count and spin.
  - create_problem — ODE problem setup wires cache, method, and model correctly.
  - evaluate_v_dot_d! — v·d NAC term computed with correct shapes and indices.
  - compute_overlap! — Overlaps remain orthonormal within tolerance.
  - calculate_Akj — Transition amplitude term Akj computed consistently with definitions.
  - evaluate_hopping_probability! — Probabilities bounded [0,1], normalized across channels.
  - execute_hop! — Electronic state and velocities updated following IESH rules.
  - algorithm comparison — Cross-integrator agreement on reference problems.
  - DecoherenceCorrectionEDC — Energy-based decoherence correction reduces overcoherence; parameters produce expected damping.

- `test/Dynamics/mdef.jl`
  - friction! — Electronic friction tensor positive semidefinite and symmetric, consistent with LDFA/EP couplings [@Gerrits2020; @Box2021; @Box2023].

- `test/Dynamics/nrpmd.jl`
  - Population correlation — NRPMD correlation consistent with normalization and symmetry [@Richardson2013; @Chowdhury2021].
  - motion! obeys Hamilton's equations — Time-derivatives satisfy Hamilton’s equations for mapping+RP variables.
  - Energy conservation $(algs[idx]) — Energy nearly conserved for symplectic integrator; Tsit5 serves as comparison baseline [@Korol2019].
  - Algorithm comparison — Agreement across integrators for observables within tolerance.
  - MInt algorithm convergence — Convergence of modified integrator order and stability as dt→0.

- `test/Dynamics/rp_ehrenfest_na.jl`
  - Algorithm comparison: $method — Consistency for ring-polymer mean-field variants on NA model.

- `test/Dynamics/rpbcme.jl`
  - BCME — Ring-polymer BCME inherits CME limits and preserves stability in bead-extended state.

- `test/Dynamics/rpcme.jl`
  - Phonon relaxation Tᵢ=$(T)T — CME relaxation rate increases with temperature and coupling strength as expected.

- `test/Dynamics/rpecmm.jl`
  - Population correlation — Mapping population autocorrelation normalized at t=0 and monotone decay.
    - Gamma = $γ — Damping γ tunes decoherence rate as designed.
  - Algorithm comparison — Reference agreement for RP mapping dynamics.

- `test/Dynamics/rpiesh.jl`
  - DynamicsVariables — RP-IESH state/caches initialized over beads.
  - create_problem — RP-IESH ODE problem constructed with beadwise couplings.
  - algorithm comparison — Consistent results across algorithms for RP-IESH test case.

- `test/Dynamics/rpmdef.jl`
  - friction! — RP friction tensor assembly correct across beads.
  - step_C! — Constrained/C step of splitting integrator updates positions/aux correctly.
  - step_O! — OU (Ornstein–Uhlenbeck) stochastic step samples noise with correct covariance [@Risken1989; @Leimkuhler2012].
  - ThermalLangevin — Thermalization drives kinetic temperature to target T.
  - MDEF — Full-step MDEF integration stable; energy/temperature behavior matches expectations.

- `test/Dynamics/spin_mapping.jl`
  - motion! obeys Hamilton's equations — Mapping-variable EOMs consistent with Hamiltonian dynamics [@He2019; @HeGong2021; @HeWu2021].
  - Energy conservation — Symplectic conservation in NVE mapping dynamics.
  - Algorithm comparison — Cross-integrator results consistent for observables.
  - MInt algorithm convergence — Convergence as dt→0 and stability over long times.

- `test/Dynamics/test_decoherence_corrections.jl`
  - apply_decoherence_correction! — Decoherence correction reduces spurious electronic coherence buildup; trends align with physical expectations.


