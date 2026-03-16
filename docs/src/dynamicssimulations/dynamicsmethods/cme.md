```@setup logging
@info "Expanding src/dynamicssimulations/dynamicsmethods/cme.md..."
start_time = time()
```
# [Surface Hopping with a Classical Master Equation (CME)](@id cme-dynamics)

Typically for surface hopping approaches such as [FSSH](fssh.md) [Tully1990](@cite) and [IESH](iesh.md) [Shenvi2009](@cite)[Gardner2023](@cite), all of the electronic states are treated explicitly. This leads to issues in IESH such as the difficulty to include relaxation of the coupled molecular-surface electronic system to equilibrium [Dou2020](@cite). An alternative approach was presented by Wenjie Dou and Joseph Subotnik in 2015 [Dou2015_II](@cite)[Dou2015_III](@cite), using the assumption that molecular electrons can be separated from the metal's electronic degrees of freedom in the weak molecule-metal interation limit [Dou2020](@cite). This allows for a Classical Master Equation to be used to describe the dynamics, where trajectories propagate classically diabatic potential energy surfaces of the molecular system, and Franck-Condon type vertical transitions at fixed nuclear positions ("hops") between surfaces account for molecule-metal interactions [Dou2015_II](@cite)[Dou2020](@cite).

The two components that dictate hopping probabilities are:
1. The hybridisation function, $\Gamma$, of the molecule-metal coupling;
2. The energy difference between the molecule (impurity) oribital and the Fermi level of the metallic system.

The Classical Master Equation propagates the probaility densities of the neutral molecules (unoccupied impurity state), $\rho_{0}$, and charged molecule (occupuied impurity state), $\rho_{1}$, within the given phase space (set of nuclear coordinates). The CME which describes the dynamics of these densities are given as:
$$
\begin{align}
\frac{\partial}{\partial t} \rho_{0} (x, p) =& \frac{\partial H_{0}(x,p)}{\partial x} \frac{\partial \rho_{0}(x,p)}{\partial p} - \frac{\partial H_{0}(x,p)}{\partial p} \frac{\partial \rho_{0}(x,p)}{\partial x} - \gamma_{0 \rightarrow 1} \rho_{0}(x,p) + \gamma_{1 \rightarrow 0} \rho_{1}(x,p) \\
\frac{\partial}{\partial t} \rho_{1} (x, p) =& \frac{\partial H_{1}(x,p)}{\partial x} \frac{\partial \rho_{1}(x,p)}{\partial p} - \frac{\partial H_{1}(x,p)}{\partial p} \frac{\partial \rho_{1}(x,p)}{\partial x} + \gamma_{0 \rightarrow 1} \rho_{0}(x,p) - \gamma_{1 \rightarrow 0} \rho_{1}(x,p)
\end{align}
$$
Where $\gamma_{0 \rightarrow 1}$ and $\gamma_{1 \rightarrow 0}$ are the hopping rates and are defined as such:
$$
\begin{align}
\gamma_{0 \rightarrow 1} =& \frac{\Gamma}{\hbar} f(h)\\
\gamma_{1 \rightarrow 0} =& \frac{\Gamma}{\hbar} (1 - f(h))
\end{align}
$$












```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
