# Docs

# Euler fluid Hamiltonian

## Hamiltonians in 1D 

### Harmonic oscilator

1-particle harmonic oscillator Hamiltonian function

$`\begin{equation}
H_{HO}(q, v) = \frac{1}{2} m v^2 + \frac{1}{2} m \omega^2 q^2,~(1)
\end{equation}`$

N-particles harmonic oscillator (independent particles, no collisions) Hamiltonian function

$`\begin{equation}
H_{NHO}(\mathbf{q}, \mathbf{v}) = \sum_i^N \frac{1}{2} m v_i^2 + \frac{1}{2} m \omega^2 q_i^2 = \frac{1}{2} m \sum_i^N v_i^2 + \frac{1}{2} m \omega^2 \sum_i^N q_i^2 \overset{_{mean}}{=} \frac{1}{2} m N \left<v^2\right> + \frac{1}{2} m \omega^2 N \left<q^2\right> ,
\end{equation}`$

Fluid harmonic oscillator ($`\overset{^{i \rightarrow \infty}_{as~\rho(x)}}{\Rightarrow}`$) Hamiltonian functional

$`\begin{equation}
\mathcal{H}_{FHO}(\rho, v) = \int \rho(x, t) \left( \frac{1}{2} v(x)^2 + \frac{1}{2} (2\pi)^2 f^2 x^2 \right) dx,~(2)
\end{equation}`$

defining Hamiltonian density

$`\begin{equation}
\mathcal{h}_{FHO}(\rho(x, t), v(x)) = \rho(x, t) \left( \frac{1}{2} v(x)^2 + \frac{1}{2} (2\pi)^2 f^2 x^2 \right),~(3)
\end{equation}`$

where linear position $`x`$ [cm], linear density $`\rho(x)`$ [g/cm], fluid velocity $`v(x) = \frac{d x}{d t}`$ [cm/s], and frequency $`f = \frac{\omega}{2\pi}`$ [rad/s].

### Ideal compressible fluid

$`\begin{equation}
\mathcal{H}_{IF}(\rho, v, \varepsilon, t) = \int \rho(x) \left( \frac{1}{2} v(x)^2 + \varepsilon(x)\right) dx,~(4)
\end{equation}`$

defining Hamiltonian density

$`\begin{equation}
\mathcal{h}_{IF}(\rho(x), v(x), \varepsilon(x)) = \rho(x) \left( \frac{1}{2} v(x)^2 + \varepsilon(x) \right),~(5)
\end{equation}`$

## Action via Hamilton-Jacobi equation in 1D

Generating function (or action) is governed by the Hamilton Jacobi equation (requires the knowledge of Hamiltonian assumed to not explicitly depend on time $`t`$)

$`\begin{equation}
\partial_t S(q, t) + H(q, \partial_q S) = 0.~(6)
\end{equation}`$.

The time-independence of the Hamiltonian allows for a separation of veriables solution

$`\begin{equation}
S(q, t) = W(q) + f(t),
\end{equation}`$ 

where $`\alpha`$ is constant.

*It needs to be highlighted, that we assume we know the dynamics of the system, i.e. coordinate $`x(t)`$ and momentum $`m v(t)`$ are known.*

The time dependent part of action is obtained from

$`\begin{equation}
\partial_t f(t) + H(q, v) = 0.~(7)
\end{equation}`$,

which leads to $`f(t) = -\alpha t`$, where $`\alpha`$ is a constant equal to the convserved energy of the system (Hamiltonian). Then the governing equation to solve for coordinate dependent part of action $`W(q)`$ is

$`\begin{equation}
H(q, \partial_q W) = - \partial_t f(t),
\end{equation}`$

which takes an explicit form depending on the form of the Hamiltonian, where we substitute momemntum $`p = m v`$ by $`\partial_q W(q)`$.

For example substituting $`m v = \partial_q W(q)`$ in the harmonic oscilator Hamiltonian (1), we get

$`\begin{equation}
\frac{1}{2m} \left( \partial_q W(q) \right)^2 + \frac{1}{2} m \omega^2 q^2 = - \partial_t f(t) .
\end{equation}`$

In case of continuous system, we write the **action functional** $`\mathcal{S}(\rho, t)`$ evaluated on the coordinate field $`\rho`$, where the Hamilton Jacobi equation for action functional reads

$`\begin{equation}
\partial_t \mathcal{S}(\rho, t) + \mathcal{H}(\rho, \delta_\rho S(\rho, t)) = 0,~(8)
\end{equation}`$

where $`\delta_\rho S`$ is the first variation of the functional $`S(\rho, t)`$ with respect to field $`\rho`$.

Similar to Hamiltonian density $`\mathcal{h}`$ we define the action density $`\mathcal{s}(\rho(x), t) = \mathcal{w}(\rho(x)) + \mathcal{f}(x, t)`$ (separation of variables thanks to conserved energy or explicitely time-independent Hamiltonian) by the following

$`\begin{equation}
\mathcal{S}(\rho, t) = \int \mathcal{s}(\rho(x), t) dx = \int \left( \mathcal{w}(\rho(x)) + \mathcal{f}(x, t) \right) dx,~(9)
\end{equation}`$

and we write the density Hamilton-Jacobi equation

$`\begin{equation}
\partial_t \mathcal{f}(x, t) + \mathcal{h}(\rho(x, t), v(x, t), \varepsilon(x, t), t) = 0,~\forall (x, t).~(10)
\end{equation}`$

It is simple to see that $`\mathcal{S}`$ defined by (9) satisfies (8) if $`\mathcal{s}`$ solves (10) $`\forall (x, t)`$, because (8) can be written as

$`\begin{equation}
\int \left( \partial_t \mathcal{s}(\rho(x, t), t) + \mathcal{h}(\rho(x, t), v(x, t), \varepsilon(x, t), t) \right) dx = 0.
\end{equation}`$

### Practical observation

One partical dynamics described by (7) leads to the solution of particle trajectory $`x(t)`$, which also solves the equivalent dynamics system of Hamilton equations $`\dot{x} = \partial_p H,~ \dot{p} = - \partial_x H`$.

Our objective is to obtain (postprocess) the action (generating function) of the dynamics (7) assuming we know the solution $`x(t)`$ and consequently $`v(t) = \dot(x)(t)`$, hence the solution $S$ of (7) only depends on time $`t`$

$`\begin{equation}
\partial_t S(t) = - H(x(t), v(t), t).~(11)
\end{equation}`$

Similarly in the case of continuous fields, our objective is to obtain (postprocess) the action functional (generating functional) of the dynamics (10) assuming we know the solution $`\rho(x, t)`$, $`v(x, t)`$, and $`\varepsilon(x, t)`$, hence the solution $\mathcal{s}$ of (10) only depends on time $`t`$ and position $`x`$

$`\begin{equation}
\partial_t \mathcal{s}(x, t) = - \mathcal{h}(\rho(x, t), v(x, t), \varepsilon(x, t), t) .~(12)
\end{equation}`$

### Action of harmonic oscilator

One particle case corresponds to governing equation (11) with one particle Hamiltonian (1) leading to

$`\begin{equation}
\partial_t S_{HO}(t) = - \frac{1}{2} m v(x(t), t)^2 - \frac{1}{2} m \omega^2 x(t)^2
\end{equation}`$

which can be easily numerically integrated (using $`S_{HO}(t=0) = 0`$).

Continuous field case corresponds to governing equation (12) with Hamiltonian density (3)

$`\begin{equation}
\partial_t \mathcal{s}_{FHO}(x, t) = - \rho(x, t) \left( \frac{1}{2} v(x)^2 + \frac{1}{2} (2\pi)^2 f^2 x^2 \right) .~()
\end{equation}`$


### Action of ideal compressible fluid

$`\begin{equation}
\partial_t \mathcal{s}_{IF}(x, t) = - \rho(x, t) \left( \frac{1}{2} v(x, t)^2 + \varepsilon(x, t) \right)
\end{equation}`$

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
