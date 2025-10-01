# Docs

# Euler fluid Hamiltonian

## Hamiltonian in 1D 
$`\begin{equation}
H(\rho, v, \varepsilon) = \int \left( \frac{1}{2}\rho(x) v(x)^2 + \rho(x) \varepsilon(x) \right) dV
\end{equation}`$

## Hamilton-Jacobi in 1D

$`\begin{equation}
\partial_t S(\rho, t) + H(\rho, v, \varepsilon, t) = 0
\end{equation}`$

Since $`\partial_t H(\rho, v, \varepsilon, t) = 0`$, $`S(q, t) = W(q) - Et`$, and $`\partial_\rho S = \partial_\rho W = \rho v`$

$`\begin{equation}
H(\rho, v, \varepsilon) = \int \left( \frac{1}{2\rho(x)}\left( \partial_\rho W(\rho(x) )\right)^2 + \rho(x) \varepsilon(x) \right) dV
\end{equation}`$

$`\begin{equation}
\int \left( \frac{1}{2\rho(x)}\left( \partial_\rho W(\rho(x) )\right)^2 \right) dV + \int \left(\rho(x) \varepsilon(x) \right) dV - E = 0
\end{equation}`$
