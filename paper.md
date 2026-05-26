---
authors:
  - name: Alexander Spears
    affiliation: '1,2'
    orcid: "0000-0002-8171-9118"
  - name: Ash Baldwin
    affiliation: '1,2'
    orcid: "0009-0008-6931-0070"
  - name: Henry Snowden
    affiliation: '2,3'
    orcid: "0009-0007-8700-7668"
  - name: Connor L. Box
    affiliation: '5'
    orcid: "0000-0001-7575-7161"
  - name: Matt Larkin
    affiliation: '2'
    orcid: "0009-0007-7547-0516"
  - name: Xuexun Lu
    affiliation: '2'
    orcid: "0009-0004-4916-5970"
  - name: Wojciech G. Stark
    affiliation: '6'
    orcid: "0000-0001-6279-2638"
  - name: James Gardner
    affiliation: '2'
    orcid: "0000-0003-1840-804X"
  - name: Nils Hertl
    affiliation: '2'
    orcid: "0000-0002-6298-3597"
  - name: Reinhard J. Maurer
    affiliation: '1,2,3,4'
    orcid: "0000-0002-3004-785X"

date: 2026-05-22

affiliations:
  - index: 1
    name: "University of Vienna, Faculty of Physics, Kolingasse 14-16, 1090 Vienna, Austria"
    ror: "03prydq77"
  - index: 2
    name: "Department of Chemistry, University of Warwick, Gibbet Hill Road, CV4 7AL, Coventry, UK"
    ror: "01a77tt86"
  - index: 3
    name: "Institute of Physical Chemistry, Georg-August University, Göttingen 37077, Germany"
    ror: "01y9bpm73"
  - index: 4
    name: "Max-Planck-Institute for Multidisciplinary Sciences, Göttingen 37077, Germany"
  - index: 5
    name: "Yusuf Hamied Department of Chemistry, University of Cambridge, , CB2 1EW Cambridge, UK"
    ror: "013meh722"
  - index: 6
    name: "Department of Chemistry, Imperial College London, White City Campus, W12 0BZ, London, UK"
    ror: "041kmwe10"

title: "NQCDynamics.jl (version 1.0): Nonadiabatic quantum classical molecular dynamics in Julia"
bibliography: main.bib
# # Remove this afterwards, I just use it for word counting
# format: 
#   pdf:
#     include-in-header: 
#       text: '\usepackage{lineno}\linenumbers'
---

# Summary
The simulation of dynamical processes of molecules and materials using classical molecular dynamics is an essential part of modern computational molecular and materials research.
However, many important and interesting scientific questions lie in the regime where the assumptions of classical molecular dynamics break down, and a first-principles quantum-mechanical treatment is computationally infeasible. 
A large variety of mixed quantum–classical dynamics has been developed to preserve key electronic quantum effects while achieving computational scaling that enables simulations of realistic, many-atom systems.[@gonzalez2021quant]

The Julia package `NQCDynamics.jl` provides an open source framework for the development of new methods to simulate non-adiabatic and quantum nuclear effects, and we have previously demonstrated the code for several different model systems [@gardner_nqcdynamicsjl_2022] in which these effects play a significant role. 
Here, we report version release 1.0, which includes code advancements aimed at making the code more suited for large-scale simulations, more accessible to users and developers, as well as more interoperable with commonly used machine learning and *ab-initio* packages.

# Statement of need
For molecular dynamics at interfaces[@DouWenjieSubotnikDynamicsAtMetalSurfaces; @bunermannElectronholePairExcitation2015; @krugerVibrationalInelasticityHighly2016a; @luoElectronholePairEffects2016], in the excited state[@gonzalez2021quant; @cigrangRoadmapMolecularBenchmarks2025a; @crespo-oteroRecentAdvancesPerspectives2018], and for dynamical systems driven by light-matter interactions [@barbattiSimulationExcitationSunlight2020; @brandbygeElectronicallyDrivenAdsorbate1995; @luntzFemtosecondLaserInduced2006], non-adiabatic effects that go beyond the Born-Oppenheimer approximation [@gonzalez2021quant] are essential to accurately describe the physical and chemical processes at play.
A full quantum description of coupled electron-nuclear dynamics is not always feasible, particularly for high-dimensional or strongly correlated condensed-phase systems, due to unfavourable computational scaling with the number of degrees of freedom.
In practice, quantum dynamics can often be approximated by treating the nuclear degrees of freedom classically while retaining a quantum mechanical description of the electronic subsystem.
These mixed quantum–classical dynamics (MQCD) methods preserve key electronic quantum effects while achieving computational scaling that enables simulations of realistic, many-atom systems.[@gonzalez2021quant]

A variety of mixed quantum-classical dynamics methods have been developed and continue to be improved both in terms of their accuracy and scale. Most commonly these methods are employed for molecular systems with few electronic states, but methods are also developed for the study of dynamics at metallic and semi-conducting surfaces.[@gardnerAssessingMixedQuantumClassical2023; @nelson_non-adiabatic_2020; @wang_recent_2016; @li_ab_2021]
In these cases, the presence of electronic bands significantly expands the scale of the problem and, in many cases, requires significant adaptations of the established methods. 

Probably, the most common approach for dynamics at metal surfaces is molecular dynamics with electronic friction (MDEF).[@head-gordon1995molec] In this framework, the electrons are treated implicitly as a bath, and the coupling between electrons and nuclear motion, represented by a friction force and a random force, is added to the conservative force present in standard molecular dynamics. MDEF in different flavours [@juaristi2008role; @maurer2016ab; @box2023ab] was successful in describing H atom scattering experiments from metal surfaces[@dorenkamp2018hydrogen; @hertl2022electronically; @box_room_2024], and vibrational dissipation.[@box2020determining] A popular alternative that treats electrons in an explicit fashion is Tully's fewest switches surface hopping (FSSH) method[@TullyMolecular1990] has been adapted to create the independent electron surface hopping method (IESH)[@shenvi09], enabling the study of a large number of electronic excitations that feature population transfer between adsorbate and metal electronic states. Another surface hopping method adapted to metal-molecule systems is the broadened classical master equation (BCME).[@dou2016broad] Both methods have been assessed against the hierarchical equations of motion (HEOM)[@tanimura1989time] method, which is a numerically exact approach to model open quantum system dynamics.[@preston2025]

# State of the field
Many of the previously described methods are not available in open-source software or only exist in specific implementations, hampering their broader adoption and reproducibility of results.
This is especially true for newly developed methods, which often have very limited documentation if a code implementation is publicly available. 
In addition, the efficiency of different implementations can vary significantly. 
As a result, benchmarking the performance of different MQCD methods for applications to a particular system poses many additional challenges.
For more mature MQCD methods such as Fewest-switches Surface Hopping, Ehrenfest dynamics or path-integral MD, high quality implementations such as [Newton-X](https://newtonx.org)[@NewtonX], [SHARC](https://sharc-md.org)[@SHARC] or [i-Pi](https://ipi-code.org)[@ipi-article] have been developed, highlighting the potential for greater adoption of these methods. 

`NQCDynamics.jl` was developed with the goal of providing a consistent foundation for the development of different methods, not unlike recent efforts in the electronic structure theory and scientific machine learning communities to establish reproducible benchmarks and open frameworks. [@althorpeEmergingOpportunitiesFuture2019a; @westermayrPerspectiveIntegratingMachine2021; @lejaeghere_error_2014; @lejaeghere_reproducibility_2016] 
By providing a library of established and developing mixed quantum-classical dynamics methods alongside a framework for deploying these methods on analytical potentials, machine learning (ML) models or electronic structure calculations, `NQCDynamics.jl` attempts to provide a platform from method development to production-scale calculations in application. 

# Software design 
Here, we present the first full release version (version 1.0) of `NQCDynamics`, released on the General Julia package registry. 

![Schematic representation of the molecular dynamics pipeline in `NQCDynamics` showing which packages of the ecosystem are involved during the respective processes.\label{fig:NQCStructure}](NQCDiagram_2.pdf){width=85%}

`NQCDynamics` attempts to strike a balance between easy entry for new users while retaining the necessary flexibility for new method developments.
This is achieved by maintaining a similar workflow across different dynamics methods and the use of consistent input and output structures based on the number of particles and degrees of freedom. 
Additions of new capabilities are made easier by splitting functionality into separate, well integrated packages (see \autoref{fig:NQCStructure}) providing electronic Hamiltonians and potential energy surfaces (`NQCModels.jl`), common methods to sample initial conditions (`NQCDistributions.jl`) and atomic structure representations (`NQCBase.jl`). 
Propagating coupled ordinary and partial differential equations requires the use of specialised (symplectic) integration algorithms using split-operator techniques. 
To enable the propagation of electronic and nuclear dynamics with different integration requirements, custom integration algorithms have been implemented in `NQCDynamics` since its inception. 
These build on the well-established `DifferentialEquations.jl` library [@rackauckas2017differentialequations], which provides an efficient framework for CPU- or GPU-based numerical integration. 
We plan to contribute the integration algorithms in `NQCDynamics.jl` to `DifferentialEquations.jl` in the future, since they have applications beyond our software package. 

The atomistic simulation ecosystem in the Julia language is quite young, and since the original release of the code, a number of different packages have become the de-facto standard.
The `NQCBase.jl` package now provides a translation layer between atomic structure representations used in `NQCDynamics` and `AtomsBase.jl`, as well as supporting the extended XYZ file format using `ExtXYZ.jl`. 
Using the `PythonCall.jl` library, we have updated our interfaces to the Python-based atomistic simulation ecosystem, allowing the use of calculators in the Atomic Simulation Environment (ASE)[@software-ase]. 
As a result, a large range of popular MLIPs such as MACE,[@batatia2022mace,Batatia2022Design] as well as electronic structure codes such as FHI-aims[@blum_ab_2009] can be used to supply the necessary potential energy surfaces for dynamics simulations.

With the version 1.0 release, NQCDynamics implements the following MQCD methods:

- classical molecular dynamics
- molecular dynamics with electronic friction[@head-gordon1995molec] (supporting multiple thermostats)
- fewest-switches surface hopping [@TullyMolecular1990]
- independent electron surface hopping (IESH) [@shenvi09]
- broadened chemical master equation dynamics [@Dou_BCME_2016]
- Ehrenfest dynamics 
- ring-polymer molecular dynamics (RPMD) [@craig2004quant]
- nonadiabatic RPMD [@chowdhury2021non-a]
- centroid ring-polymer surface hopping (RPSH) [@shushkov2012ring]
- centroid ring-polymer Ehrenfest dynamics
- generalised spin mapping approach
- extended classical mapping model [@he2021negat]

# Research impact statement

: Examples for the use of NQCDynamics.jl in various applications, including the respective contributions by each work to development. []{label="tab:nqcresearch"}

| Reference | Dynamics method used | Development contribution to NQCDynamics v1.0 |
| --- | ---- | ---- |
| [@stark_nonadiabatic_2025] | molecular dynamics with electronic friction | Interfaces to Python- and Julia-based machine learning models |
| [@box_room_2024] | molecular dynamics with electronic friction | Post-processing tools for trajectory classification |
| [@spears2026rolea] | molecular dynamics with electronic friction | Support for multiple thermostats |
| [@gardnerEfficiente2023] | IESH | Method implementation, performance improvements | 
| [@lu_H/Ge_2025] | IESH, Ehrenfest | Discretisation schemes for Anderson-Haldane models | 
| [@gardnerAssessingMixedQuantumClassical2023] | multiple methods | Benchmarking of different methods on analytical model systems |  

Since its initial release, the package has been used extensively to model nonadiabatic systems, as shown in table \ref{tab:nqcresearch}. 
A variety of MQCD methods in `NQCDynamics.jl` have been benchmarked on a range of analytical model systems.[@gardnerAssessingMixedQuantumClassical2023; @gardnerEfficiente2023]
In particular, our implementation of IESH achieves a nominal scaling of $N^3$ where $N$ is the number of electronic states, as shown in \autoref{fig:ieshscaling}. 

![Strong scaling test of a single IESH simulation trajectory with varying number of bath states. All calculations were performed on a single core of an AMD EPYC$^\text{TM}$ 7742 processors with 256GB of memory available.\label{fig:ieshscaling}](IESH_performance_plot.pdf){width=66%}

The [documentation page](https://nqcd.github.io/NQCDynamics.jl/dev/) for `NQCDynamics` and related packages has been rewritten and extended to lower the barrier to entry for new users and method developers. 
This includes automatically generated documentation based on code comments as well as dedicated pages on, e.g. setting up initial conditions or use of a certain dynamics method. 
The included tutorial examples have been updated to include newer and more relevant models that have been published since the original release of the code. 
Furthermore, a description of complete dynamics workflows used in published works by our group as well as other in-depth applications of our code is now available in our GitHub repository [NQCRecipes](https://nqcd.github.io/NQCRecipes/).

# AI Usage
Large Language Models (LLMs) were used for code development (e.g. GitHub Copilot, Claude Code) and drafting documentation changes. All code changes and documentation changes were reviewed by the authors. 
Similarly, LLMs were used to assist in the generation of parts of the manuscript, however all text was reviewed and edited by the authors.

# Acknowledgements

The authors acknowledge support via a UKRI Future Leaders Fellowship \[MR/X023109/1\], a UKRI Frontier grant \[EP/X014088/1\], an MSCA postdoctoral fellowship \[EP/Z001498/1\], a Leverhulme Trust Research Project grant \[RPG-2019-078\], and an Alexander-von-Humboldt Professorship. X.L. and M.L. are supported by EPSRC Doctoral Training Partnership studentships.

Computational resources were provided by the EPSRC-funded HEC Materials Chemistry \[EP/L000202/1, EP/R029431/1\] and the HPC-CONEXS \[EP/X035514/1\] consortia for access to the ARCHER2 UK National Computing Service, and the  
EPSRC-funded HPC Midlands+ computing centre for access to Sulis \[EP/P020232/1\].

# References
