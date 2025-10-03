# Docs

# Euler fluid Hamiltonian

## Hamiltonians in 1D 

### Harmonic oscilator

1-particle harmonic oscillator Hamiltonian function

$`\begin{equation}
H_{HO}(x, v) = \frac{1}{2} m v^2 + \frac{1}{2} m \omega^2 x^2,
\end{equation}`$

N-particles harmonic oscillator (independent particles, no collisions) Hamiltonian function

$`\begin{equation}
H_{NHO}(x, v) = \sum_i^N \frac{1}{2} m v_i^2 + \frac{1}{2} m \omega^2 x_i^2 = \frac{1}{2} m \sum_i^N v_i^2 + \frac{1}{2} m \omega^2 \sum_i^N x_i^2 = \frac{1}{2} m N \left<v^2\right> + \frac{1}{2} m \omega^2 N \left<x^2\right> ,
\end{equation}`$

Fluid harmonic oscillator ($`\overset{^{i \rightarrow \infty}_{as~\rho(x)}}{\Rightarrow}`$) Hamiltonian functional

$`\begin{equation}
H_{HFO}(\rho(x), v(x)) = \int \rho(x) \left( \frac{1}{2} v(x)^2 + \frac{1}{2} (2\pi)^2 f^2 x^2 \right) dx,
\end{equation}`$

where linear position $`x`$ [cm], linear density $`\rho(x)`$ [g/cm], fluid velocity $`v(x) = \frac{d x}{d t}`$ [cm/s], and frequency $`f = \frac{\omega}{2\pi}`$ [rad/s].

### Ideal compressible fluid

$`\begin{equation}
H_{IF}(\rho, v, \varepsilon) = \int \rho(x) \left( \frac{1}{2} v(x)^2 + \varepsilon(x) \right) dx
\end{equation}`$

## Hamilton-Jacobi



We assume that $`\partial_{\nabla \rho} W = 0`$, we get variational derivative $`\delta_\rho W = \partial_\rho W`$.


$`\begin{equation}
\partial_t S(\rho, t) + H(\rho, v, \varepsilon, t) = 0
\end{equation}`$

where Hamiltonian is defined through the hamiltonian density $`h(x)`$ (local quantity)

https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/16%3A_Analytical_Formulations_for_Continuous_Systems/16.04%3A_The_Hamiltonian_density_formulation_for_continuous_systems

$`\begin{equation}
H(\rho, v, \varepsilon, t) = \int h(\rho, v, \varepsilon, t) dV.
\end{equation}`$

For example the ideal fluid's Hamiltonian density (corresponding to total energy density) reads

$`\begin{equation}
h(x, t) = h(rho(x, t), v(x, t), \varepsilon(x, t)) = \rho(x, t) \left( v(x, t)^2 + \varepsilon(x, t)\right).
\end{equation}`$

Can we write

$`\begin{equation}
\int \left( \partial_t s(\rho, t) + h(\rho, v, \varepsilon, t) \right) dV = 0 \Leftarrow  \partial_t s(x, t) + h(x, t)
\end{equation} = 0,~ s(x, 0) = ?`$

and evaluate action functional based on the action density $`s(x)`$ (local quantity) as

$`\begin{equation}
S(\rho, t) = \int s(\rho(x, t), t) dV,
\end{equation}`$

?

Since $`\partial_t H(\rho, v, \varepsilon, t) = 0`$, $`S(q, t) = W(q) - Et`$, and $`(\delta_\rho S)(x) = (\partial_\rho S)(x) = (\partial_\rho W)(x) = \rho(x) v(x)`$

$`\begin{equation}
H(\rho, v, \varepsilon) = \int \left( \frac{1}{2\rho(x)}\left( \nabla \partial_\rho W(\rho(x) )\right)^2 + \rho(x) \varepsilon(x) \right) dV
\end{equation}`$

$`\begin{equation}
\int \left( \frac{1}{2\rho(x)}\left( \nabla \partial_\rho W(\rho(x))\right)^2 \right) dV + \int \left(\rho(x) \varepsilon(x) \right) dV - E = 0
\end{equation}`$

$`\begin{equation}
\int \left( \frac{1}{2\rho(x)}\left( \nabla \partial_\rho W(\rho(x) )\right)^2 + \rho(x) \varepsilon(x) - E \right) dV = 0
\Rightarrow \frac{1}{2\rho(x)}\left( \nabla \partial_\rho W(\rho(x) )\right)^2 + \rho(x) \varepsilon(x) - E = 0
\end{equation}`$

$`\begin{equation}
\nabla \partial_\rho W(\rho(x) ) = \sqrt{2 \rho(x) E  - 2 \rho(x)^2 \varepsilon(x)}
\end{equation}`$

$`\begin{equation}
\nabla \partial_\rho W = \sqrt{2 \rho E  - 2 \rho^2 \varepsilon}
\end{equation}`$
