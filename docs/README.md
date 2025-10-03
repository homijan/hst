# Docs

# Euler fluid Hamiltonian

## Hamiltonian in 1D 

Harmonic oscillator

$`\begin{equation}
H_{HO}(p, q) = \frac{p^2}{2 m} + \frac{1}{2} m \omega^2 q^2,
\end{equation}`$

Ideal fluid

$`\begin{equation}
H_{IF}(\rho, v, \varepsilon) = \int \left( \frac{1}{2}\rho(x) v(x)^2 + \rho(x) \varepsilon(x) \right) dV
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
