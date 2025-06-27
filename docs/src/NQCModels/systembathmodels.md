# System-bath models
This documentation is a collaborative work jointly authored by [Xuexun Lu (Hok-seon)](https://louhokseson.github.io), [Matt Larkin](https://uk.linkedin.com/in/matthew-larkin-773307219) and [Nils Hertl](https://scholar.google.com/citations?user=I65S8ZQAAAAJ&hl=en).





### Introduction

When dealing with molecule-surface systems, wherein the surface is a metal or semi-conductor; the surface can be considered as an environmental bath that allows the dissipation of energy through coupling between the molecular "system" and the bath.

The effect that the bath has on the system is described by the bath spectral density, $J(\varepsilon)$, which can take various forms dependent on the method of energy coupling.

<!-- >This is specifically for Newns-Anderson (and WideBandBath), Spin Boson samples frequency and state coupling values from the discretised bath instead of electronic energies.
>$$
>J(\varepsilon) = \int_{a}^{b} d\varepsilon' \left| V(\varepsilon') \right|^{2} \delta(\varepsilon - \varepsilon') = \left| V(\varepsilon) \right|^{2}
>$$
>Where $V(\varepsilon)$ is the "coupling function" which describes how the bath and system should be related at a given energy. -->



## Spin-Boson model
The Spin-Boson model consists of a 2 state system coupled to a bath of harmonic oscillators that treat the energy dissipation. The potential of the system is given below.
$$
V(\mathbf{\hat{R}}) = 
\begin{pmatrix}
\varepsilon + \mathbf{c}^{T} \mathbf{\hat{R}} & \Delta\\
\Delta & - \varepsilon - \mathbf{c}^{T} \mathbf{\hat{R}}
\end{pmatrix} + \frac{1}{2} \mathbf{\hat{R}}^{T} \mathbf{\Omega}^{2} \mathbf{\hat{R}} 
$$
The set of coupling coefficients ($\mathbf{c} = \{c\}$) and frequencies ($\mathbf{\Omega} = \{\omega\}$) of the bath harmonic modes required for this model are sampled from a discretisation of the bath spectral density, $J(\omega)$.

When defining the Spin-Boson model, 4 arguments are provided to the function:
```julia
using NQCModels
SpinBoson(density::SpectralDensity, N::Integer, ϵ, Δ)
```
The first 2 arguments given detail the bath discretiation, the first being the type of bath, the second (`N`) indicates the number of discretised bath modes. From these arguements, the bath discretisation functions can return the set of frequencies , $\Omega$ and coupling coefficients for each harmonic mode $j$.

The final 2 arguments `ϵ` and `Δ` define information about the 2 state system coupling to the bath. `ϵ` provides the energy bias between the two states and `Δ` the coupling between them (also referred to as the tunneling matrix element) [3].

Discretisations of the spectral density function have be implemented for two of the bath types that are prevelant in literature.
 - **Ohmic** bath (`OhmicSpectralDensity()`)
 - **Debye** bath (`DebyeSpectralDensity()`, `AltDebyeSpectralDensity`)

These functions are containers that store the relevant information needed to build a discretised bath. The actual construction and subsequent sampling is only done when `SpinBoson()` is called.

The following section details each bath discretisation method with examples of their implementation.

### Ohmic bath

The Ohmic bath spectral density function (for $\omega \geq 0$) takes the form:
$$
J(\omega) = \frac{\pi}{2} \alpha \omega e^{-\omega / \omega_{c}}  
$$
Where $\omega_{c}$ is the characteristic frequency of the bath and $\alpha$ is the Kondo parameter. These are provided as inputs to the function `OhmicSpectralDensity(ωᶜ,α)`.

The set of frequencies, $\mathbf{\Omega}$, and coupling coefficients, $\mathbf{c}$, are generated for the number of discretised bath modes, $N_{b}$ provided.
$$
\begin{align}
   \omega_{j} &= - \omega_{c} \ln\left[ 1 - \frac{j}{(1 + N_{b})} \right] \\
   c_{j} &= \sqrt{\frac{\alpha \omega_{c}}{N_{b} + 1}} \omega_{j}
\end{align}
$$

Where $j=1,...,N_{b}$.

### Debye bath
The Debye bath spectral density function takes the form:
$$
J(\omega) = 2 \lambda \frac{\omega_{c} \omega}{\omega_{c}^{2} + \omega^{2}} 
$$
Where $\omega_{c}$ is the characteristic frequency of the bath and $\lambda$ is the reorganisation. These are provided as inputs to the function `DebyeSpectralDensity(ωᶜ,λ)`.

The set of frequencies, $\mathbf{\Omega}$, and coupling coefficients, $\mathbf{c}$, are generated for the number of discretised bath modes, $N_{b}$ provided.
$$
\begin{align}
   \omega_{j} &= \omega_{c} \tan\left( \frac{\pi}{2} \left(1 - \frac{j}{(1 + N_{b})} \right) \right) \\
   c_{j} &= \sqrt{\frac{2 \lambda}{N_{b} + 1}} \omega_{j}
\end{align}
$$

Where $j=1,...,N_{b}$.

<!-- This fucntion is in the code but hasn't been exported. -->
<!-- An alternative Debye bath discretisation has also been implemented that includes an additional cuttoff frequency, `ωᵐ`. The inputs to this function are `AltDebyeSpectralDensity(ωᶜ,λ,ωᵐ)`.

The set of frequencies and coupling coefficients are again generated for the number of discretised bath modes, $N_{b}$ provided.
$$
\begin{align}
   \omega_{j} &= \omega_{c} \tan\left(j \frac{\arctan\left(\omega_{m}/\omega_{c}\right)}{N_{b}} \right) \\
   c_{j} &= \sqrt{\frac{4 \lambda \arctan\left(\omega_{m}/\omega_{c}\right)}{\pi N_{b}}} \omega_{j}
\end{align}
$$

Where $j=1,...,N_{b}$. -->

## Newns-Anderson model
The **Anderson Impurity Model (AIM)** describes a localized impurity state interacting with a continuous band of bath states. AIM is a foundamental model in condensed matter physics and quantum chemistry, introduced by P.W. Anderson in 1961. The **Newns-Anderson model** is a generalization of the AIM, which includes the possibility of **multiple impurity states** and a more complex interaction with the bath. A key advantage of using the AIM lies in its ability to yield **analytical solutions** for the energy level distribution and the hybridization (coupling) density, making it a powerful tool for theoretical analysis.
#### Theory
The **total Hamiltonian** $H$ of the AIM can be written as:
$$
H = H_{S} + H_{B} + H_{C}
$$
##### Impurity

The **Impurity Hamiltonian** describes the localized orbitals:

$$H_{S} = \sum_{i} \varepsilon_{i} d_{i}^{\dagger} d_{i}$$

Here:
* $\varepsilon_{i}$ is the energy of the $i$-th impurity orbital.
* $d_i^{\dagger}, d_i$ : creation/annihilation operators for orbital i

For the specific case of a **two-state impurity system**, the system Hamiltonian can be expressed as:

$$H_{S} = h \cdot d^{\dagger} d + U_0$$

In this two-state representation:
* $U_0$: energy of state 0
* $U_1$: energy of state 1
* $h$: the energy gap $U_1 -U_0$

##### Bath 
Describes a reservoir of states which are non-interacting and intrisically existing. The bath states can be phonons or electrons in surfaces. Usually written as
$$
H_B = \sum_{k} \epsilon_k c_k^{\dagger} c_k
$$

- $\epsilon_k$: energy of bath state $k$  
- $c_{k}^\dagger, c_{k}$: creation/annihilation operators for state $k$

##### Interaction

Captures the coupling between the impurity and the bath states. The interaction Hamiltonian can be expressed as:
$$
H_{C} = \sum_{i,k} V_{ik} d_i^{\dagger} c_k + V_{ik}^{*} c_k^{\dagger} d_i
$$

where $V$ stands for the coupling strength.

### Examples
To build a Newns-Anderson model in NQCModels.jl, you can use the [`AndersonHolstein`](@ref):
```julia
using NQCModels
...
NAmodel = AndersonHolstein(diabaticmodel, bath; couplings_rescale)
```
The `diabaticmodel` belongs to [`DiabaticModels`](@ref) is a system model describling the system Hamiltonian $H_S$ and the `bath` is a set of energy states (and their couplings to the system) that represent the environmental bath, generated by some discretisation scheme. The [`couplings_rescale`](@ref AndersonHolstein) is a scalar parameter that rescales the coupling strengths to the bath states. And the `bath` is a collection of bath states.
### Spin-Boson model
PS: Hokseon not sure whether he is the best person to write this section


## Discretisation of bath 
Some mixed quanutm classical dynamics methods that deal with system-bath simulations require the discretisation of the bath spectral density, $J(\varepsilon)$, into a finite number of discrete energy levels such that the individual state couplings are excplicitly considered during propagation. i.e. [Independent electron surface hopping](@ref adiabaticIESH) [1] [2].

The bath spectral density here is given by the following integral:
$$
J(\varepsilon) = \int_{a}^{b} d\varepsilon' \left| V(\varepsilon') \right|^{2} \delta(\varepsilon - \varepsilon') = \left| V(\varepsilon) \right|^{2}
$$
Where $V(\varepsilon)$ is the "coupling function" which describes how the bath and system should be related at a given energy [4].
Discretising this over a set of $N$ bath states gives the following form:
$$
J^{discr}(\varepsilon) = \sum_{n=1}^{N} \left|V_{n}\right|^{2} \delta \left( \varepsilon - \varepsilon_{n} \right)
$$
Where $V_{n}$ is the coupling contribution at state n. 
 
#### Discretisation under wide band limit
The **Wide Band Limit (WBL)** is a common approximation which considers the bath spectral density to be constant over a wide range of energies.
This allows in the case of direct discretisation of the bath spectral density function, $J(\varepsilon)$, that the coupling associated with each state is simply equivalent to the energy spacing between itself and its neighbouring state.

#### Discretisation with constant spacing
The simpliest choice for discretising the bath spectral density is to represent the bath as a set of `N` evenly spaced energy states according to a trapezoidal integration rule. The discretised coupling elements, $V_{n}$, are therefore given by:
$$
\left| V_{n} \right|^{2} = \left| V(\varepsilon_{n}) \right|^{2} \Delta \varepsilon_{n}
$$
Where:
 - $\Delta \varepsilon_{n}$ is the spacing between the $n$ discretised energy states 
 - $V(\varepsilon_{n})$ is the coupling function evaluated at the $n$ discretised energy states
 
When making the Wide Band Limit approximation however, the coupling function is independent of energy, $\varepsilon$. Taking $V(\varepsilon_{n})$ as an arbitrary constant allow the simplification:
$$
\left| V_{n} \right|^{2} \approx \Delta \varepsilon_{n}
$$

This is implemented as `TrapezoidalRule()` which takes the following arguements:
```julia
using NQCModels
TrapezoidalRule(M, bandmin, bandmax)
```
Where:
 - `M` = number of discretised bath states to be generated
 - `bandmin` / `bandmax` = miniumum / maximum energies that define the bath energy range

The discretised bath energy states, $\varepsilon_{n}$, and coupling elements, $V_{n}$, are calculated as follows:
$$
\begin{align*}
   \varepsilon_{n} &= a + \frac{(n - 1)(b - 1)}{M} \\
   V_{n} &= V(\varepsilon_{n})\sqrt{(b-a)/M} = \Delta \varepsilon_{n}
\end{align*}
$$ 
Where:
 - $n$ = index of the discretised state
 - $a$ = minimum energy of bath discretisation (`bandmin`)
 - $b$ = maximum energy of bath discretisation (`bandmax`)

This discretisation method is effective and given enough states will always accurately represent the spectral density, but methods such as IESH using this discretisation scheme will suffer from long computation times as the number of bath states $N$ increases. 

>[!note]
Typically keeping $N<100$ for this discretisation scheme is sensible.

#### Discretisation by Gauss-Legendre quadrature


#### Discretisation with a gap
Two numerical discretisation methods which performs a band gap in the middle of the continuum are introduced and named as [`GapTrapezoidalRule`](@ref) and [`GapGaussLegendre`](@ref).


##### Gapped Trapezoidal Rule
The *Gapped Trapezoidal* Rule discretises a band into two evenly spaced continuums, which is particularly useful for systems with a band gap, such as semiconductors.

The energy $\epsilon_k$ for each discretised state $k$ is defined as:

$$ \epsilon_k = \begin{cases}E_{\text{F}} - \Delta E + (k-1) \times \frac{\Delta E - E_{\text{gap}}}{M} &\text { if } k \leq M/2 \\ E_{\text{F}} + \frac{E_{\text{gap}}}{2} + (k-M/2 -1) \times \frac{\Delta E - E_{\text{gap}}}{M}&\text { otherwise }\end{cases} $$

These states are associated with constant coupling weights $\omega_k$:

$$\omega_k = \frac{\Delta E - E_{\text{gap}}}{M}$$

where:
* $\Delta E$ is the total energy range of the bath.
* $E_{\text{gap}}$ is the size of the gap located in the middle of the band.
* $E_{\text{F}}$ is the Fermi energy.
* $M$ is the total number of discretised states.
##### Gapped Gauss Legendre 
The *Gapped Gauss Legendre* discretisation method provides a robust way to sample the bath spectral density, specifically designed to incorporate a band gap in the middle of the continuum. This approach leverages Gauss-Legendre quadrature for accurate representation.

The discretisation for the energy $\epsilon_k$ of each state $k$ is defined as follows:

$$
\epsilon_k = \begin{cases}E_{\text{F}} - \left[\frac{\Delta E-E_{\text{gap}}}{2} (1+x_{\text L,M/2-k+1}) + E_{\text{gap}}\right]/2 &\text { if } k \leq M/2 \\
E_{\text{F}} + \left[\frac{\Delta E-E_{\text{gap}}}{2} (1+x_{\text L,k-M/2}) + E_{\text{gap}}\right]/2&\text { otherwise }\end{cases}
$$

with conjugate weights:

$$
w_k = \begin{cases}
\Delta E w_{\text{L},M/2-k+1} /2 &\text { if } k \leq M/2 \\
\Delta E w_{\text{L},k-M/2} &\text { otherwise }
\end{cases}
$$

**Parameters:**

* $x_{\text{L},i}$ and $w_{\text{L},i}$: These represent the knot points and weights, respectively, obtained from Legendre quadrature over the interval $[-1,1]$. For this formulation, $M/2$ knots are used.
* $E_{\text{F}}$: The Fermi level, which is typically set to 0 eV or the centre of the continuum band.
* $M$: The total number of discretisation points. It is crucial that $M$ is an even number for this formulation to be valid.
* $\Delta E$: The total width of the continuum band.
* $E_{\text{gap}}$: The size of the band gap.



##### Gapped discretisation examples
To build a discretised bath with a gap in the middle, you can use the [`GapTrapezoidalRule`](@ref) or [`GapGaussLegendre`](@ref) methods.
<p align="center">
  <img src="../assets/system-bath-model/DOS_bath_discretisation_compare_test.svg" alt="Discretisation with a gap" style="display: block; margin: 0 auto;">
</p>

The Julia script to reproduce the above figure is available in the [plot_bath_DOS.jl](../assets/system-bath-model/plot_bath_DOS.jl).


#### Windowed Discretisation
> - `WindowedTrapezoidalRule()` implementation
> - Detail WTR8 with coupling correction removed
> - This method allows for known energy regions of interest to be finely sampled without the detriment of long computation times that would arise from increasing the number of states in the discretisation
> - Benefit when compared to a Gauss-Legendre discretisation scheme is that you can finely discretise a larger energy window than what is allowed by a Gauss-Legendre discretisation
> - This discretisation method is important if you want to sample your electronic state occupations from a non-equilibrium electorn distribution, where there is more than one region of interest in the probability distribution (not just the Fermi Level).  

#### How to implement and develop your own discretisation


## References

[1] N. Shenvi, S. Roy, and J. C. Tully, "Nonadiabatic dynamics at metal surfaces: Independent-electron surface hopping", *J. Chem. Phys.* **130**, 174707 (2009). [DOI: 10.1063/1.3125436](https://doi.org/10.1063/1.3125436)

[2] J. Gardner, D. Corken, S. M. Janke, S. Habershon, and R. J. Maurer, "Efficient implementation and performance analysis of the independent electron surface hopping method for dynamics at metal surfaces", *J. Chem. Phys.* **158**, 064101 (2023). [DOI: 10.1063/5.0137137](https://doi.org/10.1063/5.0137137)

[3] X. He and J. Liu, "A new perspective for nonadiabatic dynamics with phase space mapping models", *J. Chem. Phys.* **151**, 024105 (2019). [DOI: 10.1063/1.5108736](https://doi.org/10.1063/1.5108736)

[4] I. de Vega, U. Schollwöck and F. A. Wolf, "How to discretize a quantum bath for real-time evolution", *Phys. Rev. B* **92**, 155126 (2015). [DOI: 10.1103/PhysRevB.92.155126](https://journals.aps.org/prb/abstract/10.1103/P) 