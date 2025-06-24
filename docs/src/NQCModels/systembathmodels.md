# System-bath models




### Introduction
PS: Hoskeon wanna say something in this section Just something enlightening for the reader to understanding the genereal idea of system bath coupling.

When dealing with molecule-surface systems, wherein the surface is a metal or semi-conductor; the surface can be considered as an environmental bath that allows the dissipation of energy through coupling between the molecular "system" and the bath.

The effect that the bath has on the system is described by the bath spectral density, $J(\varepsilon)$, which can take various forms dependent on the method of energy coupling.

>This is specifically for Newns-Anderson (and WideBandBath), Spin Boson samples frequency and state coupling values from the discretised bath instead of electronic energies.
>$$
>J(\varepsilon) = \int_{a}^{b} d\varepsilon' \left| V(\varepsilon') \right|^{2} \delta(\varepsilon - \varepsilon') = \left| V(\varepsilon) \right|^{2}
>$$
>Where $V(\varepsilon)$ is the "coupling function" which describes how the bath and system should be related at a given energy.



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
SpinBoson(density::SpectralDensity, N::Integer, ϵ, Δ)
```
The first 2 arguments given detail the bath discretiation, the first being the type of bath, the second (`N`) indicates the number of discretised bath modes. From these arguements, the bath discretisation functions can return the set of frequencies , $\Omega$ and coupling coefficients for each harmonic mode $j$.

The final 2 arguments `ϵ` and `Δ` define information about the 2 state system coupling to the bath. `ϵ` provides the energy bias between the two states and `Δ` the coupling between them (also referred to as the tunneling matrix element).

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

PS: Hokseon think we can briefly introduce the theory and link one or two of our papers for furtther reading for the interested reader. Bring more visitations/downloads to our papers!



### Discretisation of bath 
talk something about why we need this. good for the IESH.

#### Discretisation with a gap
Hokseon'job
#### Your dicrestisation @mlarkin863

#### How to implement and develop your own discretisation

