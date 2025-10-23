# Docs

# Re-engineering of action functional from Lagrangian density

The action functional is defined as

$`\begin{equation}
S[f](\tilde{\tau}) = \int_{0}^{\tilde{\tau}} L[f] d\tau = \int_{0}^{\tilde{\tau}} \int_\Omega \mathcal{L}(f(a, \tau), \dot{f}(a, \tau), \partial_a f(a, \tau), a, \tau) d^3a~d\tau,~(L0)
\end{equation}`$

where *Lagrangian desnity* $`\mathcal{L}`$ is a function of the field and its derivatives at given point $`(a, \tau)`$. 

The integral formulation of action $`(L0)`$ can be written in differential form as

$`\begin{equation}
\frac{d S[f]}{d \tau}(\tau) = L[f](\tau) = \int_\Omega \mathcal{L}(f(a, \tau), \dot{f}(a, \tau), \partial_a f(a, \tau), a, \tau) d^3a.~(L1)
\end{equation}`$

The *conjugate momentum* to $`f`$ is defined as

$`\begin{equation}
\pi = \partial_{\dot{f}} \mathcal{L}.~(L2)
\end{equation}`$

## Proposal of action density equivalence

The action density field $`\mathcal{S}_f(a, t)`$ solving the following ordinary differential equation for given field $`f`$

$`\begin{equation}
\frac{d \mathcal{S}_f}{d \tau}(a, \tau) = \mathcal{L}(f(a, \tau), \dot{f}(a, \tau), \partial_a f(a, \tau), a, \tau),~(L3)
\end{equation}`$

leads to the same action as evaluated by $`(L0)`$ written mathematically as

$`\begin{equation}
S[f](\tilde{\tau}) = \int_\Omega \mathcal{S}_f(a, \tilde{\tau}) d^3a.~(L4)
\end{equation}`$

## Ideal fluid example

Let's $`a`$ is the fluid coordinate, then the ideal fluid Lagrangian density function reads

$`\begin{equation}
\mathcal{L}^{if}(\rho_0(a), \dot{f}(a, \tau), \partial_a f(a, \tau)) = \rho_0(a) \left( \frac{1}{2}\dot{f}(a, \tau)^2 - \varepsilon\left( s_0(a), \frac{\rho_0(a)}{\partial_a f(a, \tau)} \right) \right),~(L5)
\end{equation}`$

where $`\varepsilon(s, \rho) = c \rho^{\gamma - 1} \exp(\alpha s)`$, where $`c`$, $`\gamma`$, and $`\alpha`$ are constants, is the internal energy potential depending on density $`\rho(a, t) = \frac{\rho_0(a)}{\partial_a f(a, \tau)}`$. Note that $`\rho_0(a)`$ (density in with respect to fluid coordinates) and $`s_0(a)`$ (isentropic processs) do not change with $`\tau`$.

The conjugate momentum of field $`f`$ defined by $`(L2)`$ from $`(L5)`$ is

$`\begin{equation}
\pi(a, \tau) = \partial_{\dot{f}} \mathcal{L}_=^{if}(a, \tau) = \rho_0(a) \dot{f}(a, \tau).~(L6) 
\end{equation}`$

Finally, we obtain the ideal fluid action density $`\mathcal{S}^{if}_f(a, \tau)`$ for a given function $`f`$ from $`(L3)`$ using ideal fluid Lagrangian density $`(L5)`$ by solving ordinary differential equation

$`\begin{equation}
\frac{d \mathcal{S}_f^{if}}{d \tau}(a, \tau) = \rho_0(a) \left( \frac{1}{2}\dot{f}(a, \tau)^2 - \varepsilon\left( \frac{\rho_0(a)}{\partial_a f(a, \tau)} \right) \right).~(L7)
\end{equation}`$

It can be shown that $`\dot{f}(a, \tau)`$ and $`\partial_a f(a, \tau)`$ are equivalent to Euler variables of ideal fluid

$`\begin{align}
\dot{f}(a, \tau) &= v(x(a, \tau), \tau),~(L8)
\\
\partial_a f(a, \tau) &= \frac{\rho_0(a)}{\rho(x(a, \tau), \tau)},~(L9)
\end{align}`$

where $`v(x, \tau)`$ and $`\rho(x, \tau)`$ are solution to Euler equations

$`\begin{align}
\rho \left( \partial_\tau v + v \partial_x v \right) &= - \partial_x p,
\\
\partial_\tau \rho + \partial_x \left( \rho v \right) &= 0,~(L10)
\\
\partial_\tau s + v \partial_x s &= 0,
\end{align}`$

with isentropic closure $`p = (\gamma - 1) \rho \varepsilon(s, \rho) = c_0 \rho^\gamma`$, $`c_0`$ constant in fluid coordinates.

Note that the realtion between Eulerian coordinate $`x`$ and Lagrangian (fluid) coordinate $`a`$ reads

$`\begin{align}
x(a, \tau) &= a + \int_0^\tau v(x(a, \tilde{\tau}), \tilde{\tau}) d\tilde{\tau},~(L11)
\\
a(x, \tau) &= x - \int_0^\tau v(x, \tilde{\tau}) d\tilde{\tau},~(L12)
\end{align}`$

$`x(a, 0)`$

# Re-engineering of action functional

Hamiltonian functional

$`\begin{equation}
\mathcal{H}[\rho, v, \varepsilon] = \int \mathcal{h}[\rho(x), v(x), \varepsilon(x)](x) dx,~(1)
\end{equation}`$

where $`\mathcal{h}(\rho(x), v(x), \varepsilon(x))`$ is the Hamiltonian density, e.g. ideal fluid
$`\mathcal{h}_{IF}(\rho(x), v(x), \varepsilon(x)) = m \rho(x) \left( \frac{1}{2} v(x)^2 + \varepsilon(x) \right)`$.

Similarly action functional

$`\begin{equation}
\mathcal{S}[\rho, t] = \int \mathcal{s}[\rho, t](x, t) dx,~(2)
\end{equation}`$

where $`\mathcal{s}(\rho(x, t), t)`$ is the action density.

The Hamilton-Jacobi equation for continuum ($`\delta_\rho \mathcal{S} \Rightarrow \partial_x \left( \partial_\rho \mathcal{s} \right)`$)

$`\begin{equation}
\mathcal{h}\left( \rho(x, t), \partial_x \left( \partial_\rho \mathcal{s}|_{\rho(x, t)} \right) \right) = - \partial_t \mathcal{s}(x, t),~(3)
\end{equation}`$

can be discretized in the case of ideal fluid Hamiltonian density as

$`\begin{equation}
\left( \frac{\frac{\mathcal{s}(x_i, t_n) - \frac{\mathcal{s}(x_i, t_{n-1}) + \mathcal{s}(x_{i-1}, t_n)}{2}}{\rho(x_i, t_n) - \frac{\rho(x_i, t_{n-1}) + \rho(x_{i-1}, t_n)}{2}} - \frac{\frac{\mathcal{s}(x_i, t_{n-1}) + \mathcal{s}(x_{i-1}, t_n)}{2} - \mathcal{s}(x_{i-1}, t_{n-1})}{\frac{\rho(x_i, t_{n-1}) + \rho(x_{i-1}, t_n)}{2} - \rho(x_{i-1}, t_{n-1})}}{x_i - x_{i-1}} \right)^2 + 2 \left( m \rho(x_i, t_n) \right)^2 \varepsilon(x_i, t_n) + 2 m \rho(x_i, t_n) \frac{\mathcal{s}(x_i, t_n) - \mathcal{s}(x_i, t_{n-1})}{\Delta t} = 0,~(4)
\end{equation}`$

allowing for **re-engineeing of action density** $`\mathcal{s}(x, t)`$.

## Lagrangian density

$`\begin{align}
S_\pi[f] &= \int_{\tau_0}^{\tau_1} L[f, \dot{f}] d\tau = \int_{\tau_0}^{\tau_1} \int_\Omega \mathcal{L}(f(a, \tau), \dot{f}(a, \tau), a, \tau) d^3a~d\tau 
\\
&\overset{^{ideal}_{fluid}}{=} \int_{\tau_0}^{\tau_1} \int_\Omega \rho_0(a) \left( \frac{1}{2}\dot{f}(a, \tau)^2 - \varepsilon(a, \tau) \right) d^3 a~d\tau \overset{_{\pi = \rho_0 \dot{f}}}{=} \int_{\tau_0}^{\tau_1} \int_\Omega \left( \frac{1}{2 \rho_0(a)} \pi(a, \tau)^2 - \rho_0(a) \varepsilon(a, \tau) \right) d^3 a~d\tau,~(L1)
\end{align}`$

where *Lagrangian desnity* $`\mathcal{L}`$ is a function of the field and its derivatives at given point $`(a, \tau)`$. Note that $`\rho_0(a)`$ does not change with $`\tau`$.

$`\delta_{f(a, \tau)} S_\pi[f] = \pi(a, \tau)`$

## Functional derivative

Given functional $`G[\rho]`$ on function $`\rho`$ defined at every point $\mathbf{s}$ of the volume domain $`\Omega`$

$`\begin{equation}
G[\rho] = \int_\Omega g(\mathbf{s}, \rho(\mathbf{s}), \nabla_{\mathbf{s}} \rho(\mathbf{s})) d\mathbf{s},~(fd1)
\end{equation}`$

the functional derivative of $`G[\rho]`$ at point $`\mathbf{s}`$ is (Formula section of https://en.wikipedia.org/wiki/Functional_derivative )

$`\begin{equation}
\delta_{\rho(\mathbf{s})} G[\rho] = \left( \partial_\rho g \right) (\mathbf{s}) - \left( \nabla_\mathbf{s} \cdot \partial_{\nabla_\mathbf{s} \rho} g \right)(\mathbf{s}).~(fd2)
\end{equation}`$

### Ideal fluid action and functional derivative

Note, that if $`\mathbf{s} = (a, \tau)`$ in $`(fd2)`$ and $`g`$ is the ideal fluid Lagrangian density from $`(L1)`$, $`\mathcal{L}(a, \dot{f}, \partial_a f) = \rho_0(a) \left( \frac{1}{2}\dot{f}^2 - \varepsilon \left( s_0(a), \frac{\rho_0(a)}{\partial_a f} \right) \right)`$, we obtain Euler-Lagrange equations (Euler–Lagrange equations section of https://en.wikipedia.org/wiki/Lagrangian_(field_theory) )

$`\begin{equation}
\partial_f \mathcal{L} = \partial_{s^i} \left( \partial_{(\partial_{s^i} f)} \mathcal{L} \right) \overset{_{\partial_f \mathcal{L} = 0}}{\Rightarrow} \partial_\tau \left( \partial_{(\partial_\tau f)} \mathcal{L} \right) + \partial_a \left( \partial_{(\partial_a f)} \mathcal{L} \right) = 0 \Rightarrow \rho_0 \partial^2_{\tau^2} f - \partial_a \left( \frac{\rho_0^2}{(\partial_a f)^2} \partial_{\frac{\rho_0(a)}{\partial_a f}} \varepsilon \right) = 0,~(L2)
\end{equation}`$

which can be found to be equivalent to second Newton law of fluid 

$`\begin{equation}
\rho \dot{v} = - \partial_x p
\end{equation}`$

where we used $`\partial_\tau = D_t = \dot{}`$ and $`v = D_t f = \partial_\tau f`$, Jacobian determinant $`J = \partial_a f`$, hence density $\rho = \frac{\rho_0}{\partial_a f}$, ideal gas equation of state $`\varepsilon = \frac{p}{(\gamma - 1) \rho}`$, hence $`\partial_{\frac{\rho_0(a)}{\partial_a f}} \varepsilon = -\frac{p}{(\gamma - 1)\rho^2}`$. Note that $`\partial_a = J \partial_x`$.

## To be digested

Considering (Goldstein3-10.12) we see that functional

$`\begin{align}
D_\tau S[f](\tau) &= \delta_{f(\tau)} S[f](\tau) D_\tau f(\tau) + \delta_\tau S[f](\tau)
\\
&\overset{_{(L1)}}{=} L[f(\tau), \dot{f}(\tau)]
\end{align}`$

## MG's notation

$`\mathcal{S}_p[f(x)](\tau)`$

$`\begin{equation}
\mathcal{H}\left[ f(x), \delta_f \mathcal{S}_p[f(x)](\tau) \right] = - \partial_\tau \mathcal{S}_p[f(x)](\tau) ,~(3)
\end{equation}`$

$`\partial_x \left( \partial_f \mathcal{S}_p[f(x)](\tau) \right)`$

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
\mathcal{H}_{FHO}(\rho, v) = \int \rho(x, t) \left( \frac{1}{2} m v(x)^2 + \frac{1}{2} m (2\pi)^2 f^2 x^2 \right) dx,~(2)
\end{equation}`$

defining Hamiltonian density

$`\begin{equation}
\mathcal{h}_{FHO}(\rho(x, t), v(x)) = \rho(x, t) \left( \frac{1}{2} m v(x)^2 + \frac{1}{2} m (2\pi)^2 f^2 x^2 \right),~(3)
\end{equation}`$

where linear position $`x`$ [cm], linear density $`\rho(x)`$ [1/cm], fluid velocity $`v(x) = \frac{d x}{d t}`$ [cm/s], and frequency $`f = \frac{\omega}{2\pi}`$ [rad/s].

### Ideal compressible fluid

$`\begin{equation}
\mathcal{H}_{IF}(\rho, v, \varepsilon, t) = \int m \rho(x) \left( \frac{1}{2} v(x)^2 + \varepsilon(x)\right) dx,~(4)
\end{equation}`$

defining Hamiltonian density

$`\begin{equation}
\mathcal{h}_{IF}(\rho(x), v(x), \varepsilon(x)) = m \rho(x) \left( \frac{1}{2} v(x)^2 + \varepsilon(x) \right),~(5)
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

Similar to Hamiltonian density $`\mathcal{h}`$ we define the action density $`\mathcal{s}(\rho(x, t), t) = \mathcal{w}(\rho(x, t)) + \mathcal{f}(x, t)`$ (separation of variables thanks to conserved energy or explicitely time-independent Hamiltonian) by the following

$`\begin{equation}
\mathcal{S}(\rho, t) = \int \mathcal{s}(\rho(x, t), t) dx = \int \left( \mathcal{w}(\rho(x, t)) + \mathcal{f}(x, t) \right) dx,~(9)
\end{equation}`$

and we write the density Hamilton-Jacobi equation, applying the variable separation,

$`\begin{equation}
\partial_t \mathcal{f}(x, t) = - \mathcal{h}(\rho(x, t), v(x, t), \varepsilon(x, t)),~(10a)
\end{equation}`$

and

$`\begin{equation}
\mathcal{h}\left( \rho(x, t), \partial_x \left( \partial_\rho \mathcal{w}|_{\rho(x, t)} \right) \right) = - \partial_t \mathcal{f}(x, t),~(10b)
\end{equation}`$

where we will use a discrete approximation of the local variation 

$`\begin{equation}
\partial_\rho \mathcal{w}|_{\rho(x, t)}
\approx
\frac{\mathcal{w}(x, t) - \mathcal{w}(x-\Delta x, t-\Delta t)}{\rho(x, t) - \rho(x-\Delta x, t-\Delta t)} .~(11)
\end{equation}`$

It is simple to see that $`\mathcal{S}`$ defined by (9) satisfies (8) if $`\mathcal{s} = \mathcal{w} + \mathcal{f}`$ solves (10) $`\forall (x, t)`$, because (8) can be written as

$`\begin{align}
0 &= \int \left( \partial_t \left( \mathcal{w}(\rho(x, t)) + \mathcal{f}(x, t) \right) + \mathcal{h}(\rho(x, t), v(x, t), \varepsilon(x, t)) \right) dx 
\\
& = \int \partial_t \left( \mathcal{w}(\rho(x, t)) + \mathcal{f}(x, t) \right) dx + \int \left( \mathcal{h}(\rho(x, t), v(x, t), \varepsilon(x, t)) \right) dx
\\
&= \partial_t \mathcal{S}(\rho, t) + \mathcal{H}(\rho, v, \varepsilon).
\end{align}`$

Equation 10b is non-obvious, but it originates from *Wasserstein gradient* required when working with continuous fields (such as density) and functional (density) Hamiltonian formulation. The formula is motivated by the following

$`\begin{multline}
\mathcal{H}\left( \rho, \delta_{\rho} \mathcal{S}(\rho, t) \right) = - \partial_t \mathcal{S}(\rho, t) \overset{^{HJD}_{Chow~(1.1)}}{\Rightarrow}
\mathcal{h}\left( \rho(x, t), \partial_x \left( \left( \left( \delta_{\rho} \mathcal{S} \right) (\rho) \right)(x) \right) \right) = - \partial_t \mathcal{f}(x, t) \overset{^{functional~derivative}_{of~integral~(9)}}{\Rightarrow} \mathcal{h}\left( \rho(x, t), \partial_x \left( \partial_\rho \mathcal{w}|_{\rho(x, t)} \right) \right) = - \partial_t \mathcal{f}(x, t)
\end{multline}`$

## Algorithm in 1D

Let's start by defining explicit, yet generally applicable, functional dependence of Hamiltonian density (slight reformulation of (5))

$`\begin{equation}
\mathcal{h}(\rho(x), v(x), \varepsilon(x)) = \frac{1}{2 m \rho(x)} \left( m \rho(x) v(x) \right)^2 + m \rho(x) \varepsilon(x).~(12)
\end{equation}`$

Hamilton-Jacobi model (10) can be rewrriten with (12) as

$`\begin{align}
\partial_t \mathcal{f}(x, t) &= - \frac{1}{2} m \rho(x, t) v(x, t)^2 - m \rho(x, t) \varepsilon(x, t),
\\
\left( \partial_x \left( \partial_\rho \mathcal{w}|_{\rho(x, t)} \right) \right)^2
&= - 2 m \rho(x, t) \left( m \rho(x, t) \varepsilon(x, t) + \partial_t \mathcal{f}(x, t) \right),
\end{align}`$

which we further simplify by substituting $`\partial_t \mathcal{f}(x, t) = - \frac{1}{2} m \rho(x, t) v(x, t)^2 - m \rho(x, t) \varepsilon(x, t)`$ into the second equation

$`\begin{align}
\partial_t \mathcal{f}(x, t) &= - \frac{1}{2} m \rho(x, t) v(x, t)^2 - m \rho(x, t) \varepsilon(x, t),~(13a)
\\
\left( \partial_x \left( \partial_\rho \mathcal{w}|_{\rho(x, t)} \right) \right)^2
&= \left( m \rho(x) v(x) \right)^2 \Rightarrow \partial_x \left( \partial_\rho \mathcal{w}|_{\rho(x, t)} \right) = m \rho(x, t) v(x, t),~(13b)
\end{align}`$

where the last equality in (13b) concludes the Hamilton-Jacobi model with unknowns $`\mathcal{f}(x, t)`$ and $`\mathcal{w}(x, t)`$.

### Discrete algorithm in 1D

Let's define a discrete approximation to $`\partial_x \left( \partial_\rho \mathcal{f}|_{\rho(x, t)} \right)`$. Along with the functional approximation (11), we define the following approximation

$`\begin{equation}
\partial_x \left( \partial_\rho \mathcal{f}|_{\rho(x_i, t_i)} \right) \approx
\frac{\frac{\mathcal{f}(x_i, t_n) - \frac{\mathcal{f}(x_i, t_{n-1}) + \mathcal{f}(x_{i-1}, t_n)}{2}}{\rho(x_i, t_n) - \frac{\rho(x_i, t_{n-1}) + \rho(x_{i-1}, t_n)}{2}} - \frac{\frac{\mathcal{f}(x_i, t_{n-1}) + \mathcal{f}(x_{i-1}, t_n)}{2} - \mathcal{f}(x_{i-1}, t_{n-1})}{\frac{\rho(x_i, t_{n-1}) + \rho(x_{i-1}, t_n)}{2} - \rho(x_{i-1}, t_{n-1})}}{x_i - x_{i-1}},~(14)
\end{equation}`$

where $`\frac{g_{in-1} + g_{i-1n}}{2}`$ geometrically coresponds to the center of $`(i, i-1) \times (n, n-1)`$ coordinate square (targeting symmetry leading to symplectic discretization).

Discrete version of (13) using the approximation (14) reads

$`\begin{align}
\frac{\mathcal{f}(x_i, t_n) - \mathcal{f}(x_i, t_{n-1})}{\Delta t} &= - \frac{1}{2} m \rho(x_i, t_n) v(x_i, t_n)^2 - m \rho(x_i, t_n) \varepsilon(x_i, t_n),
\\
\frac{\frac{\mathcal{w}(x_i, t_n) - \frac{\mathcal{w}(x_i, t_{n-1}) + \mathcal{w}(x_{i-1}, t_n)}{2}}{\rho(x_i, t_n) - \frac{\rho(x_i, t_{n-1}) + \rho(x_{i-1}, t_n)}{2}} - \frac{\frac{\mathcal{w}(x_i, t_{n-1}) + \mathcal{w}(x_{i-1}, t_n)}{2} - \mathcal{w}(x_{i-1}, t_{n-1})}{\frac{\rho(x_i, t_{n-1}) + \rho(x_{i-1}, t_n)}{2} - \rho(x_{i-1}, t_{n-1})}}{x_i - x_{i-1}} &= m \rho(x_i, t_n) v(x_i, t_n),
\end{align}`$

which can be formulated in its forward advection form

$`\begin{align}
\mathcal{f}(x_i, t_n) &= \mathcal{f}(x_i, t_{n-1}) - \Delta t \left( \frac{1}{2} m \rho(x_i, t_n) v(x_i, t_n)^2 +  m \rho(x_i, t_n) \varepsilon(x_i, t_n) \right),~(d13a)
\\
\mathcal{w}(x_i, t_n) &= \frac{\mathcal{w}(x_i, t_{n-1}) + \mathcal{w}(x_{i-1}, t_n)}{2} + \left( \rho(x_i, t_n) - \frac{\rho(x_i, t_{n-1}) + \rho(x_{i-1}, t_n)}{2} \right) \left( \frac{\frac{\mathcal{w}(x_i, t_{n-1}) + \mathcal{w}(x_{i-1}, t_n)}{2} - \mathcal{w}(x_{i-1}, t_{n-1})}{\frac{\rho(x_i, t_{n-1}) + \rho(x_{i-1}, t_n)}{2} - \rho(x_{i-1}, t_{n-1})} + (x_i - x_{i-1}) m \rho(x_i, t_n) v(x_i, t_n) \right).~(d13b)
\end{align}`$

Note that the action density model (10) uses separation of variables if $`\partial t \mathcal{H} = 0`$. This contraint can be relaxed and (10) would take the common form of Hamilton-Jacobi equation for action density

$`\begin{equation}
\mathcal{h}\left( \rho(x, t), \partial_x \left( \partial_\rho \mathcal{s}|_{\rho(x, t)} \right) \right) = - \partial_t \mathcal{s}(x, t),~(15)
\end{equation}`$

which can be discretized as

$`\begin{equation}
\left( \frac{\frac{\mathcal{s}(x_i, t_n) - \frac{\mathcal{s}(x_i, t_{n-1}) + \mathcal{s}(x_{i-1}, t_n)}{2}}{\rho(x_i, t_n) - \frac{\rho(x_i, t_{n-1}) + \rho(x_{i-1}, t_n)}{2}} - \frac{\frac{\mathcal{s}(x_i, t_{n-1}) + \mathcal{s}(x_{i-1}, t_n)}{2} - \mathcal{s}(x_{i-1}, t_{n-1})}{\frac{\rho(x_i, t_{n-1}) + \rho(x_{i-1}, t_n)}{2} - \rho(x_{i-1}, t_{n-1})}}{x_i - x_{i-1}} \right)^2 + 2 \left( m \rho(x_i, t_n) \right)^2 \varepsilon(x_i, t_n) + 2 m \rho(x_i, t_n) \frac{\mathcal{s}(x_i, t_n) - \mathcal{s}(x_i, t_{n-1})}{\Delta t} = 0,~(16)
\end{equation}`$

which is a more general version (though nonlinear) of (d13), because it allows time dependent Hamiltonian.

Note, that $`\varepsilon(x_i) = \frac{1}{2} \omega^2 x_i^2`$ in the case of *continous* harmonic oscilator.

**Exercise: Derive action of harmonic oscilator from its known dynamics**

The Hamilton-Jacobi equation (6) for harmonics oscilator via its Hamitlonian (1) reads

$`\begin{equation}
\partial_t S(q, t) + \frac{1}{2 m} \left( \partial_q S \right)^2 + \frac{1}{2} m \omega^2 q^2 = 0,~(17)
\end{equation}`$

which can be rewritten using separation of variables $`S(q, t) = W(q) + f(t)`$ as

$`\begin{equation}
\left( \partial_q W \right)^2 = - m^2 \omega^2 q^2 - 2 m \partial_t f(t) \overset{_{\partial_t f(t) = H(q, v)}}{=} (m v)^2 \Rightarrow \partial_q W = mv.~(18)
\end{equation}`$

The discrete version of the last equality of (18) gives

$`\begin{equation}
\frac{W_i - W_{i-1}}{q_i - q_{i-1}} = m v_i,~(19)
\end{equation}`$

where $`W_i`$ and $`v_i`$ corresponds to value at $`q_i`$. Further, the discrete coodrinate $`q_i = q(t_i)`$ correspond to discrete time $`t_i`$ based on the known dynamics of harmonic oscillator. 

The equivalent *continous* formula for the fluid harmonic oscilator (d13b) can be reduced considering $`\rho(x_i, t_n) = \delta_{in}`$, where $`i`$ is the spatial index and $`n`$ is the temporal index, hence

$`\begin{equation}
\mathcal{w}(x_i, t_n) = \mathcal{w}(x_{i-1}, t_{n-1}) + ( x_i - x_{i-1}) m v(x_i, t_n),~(20)
\end{equation}`$

which concludes the motivation of this excerices that (20) is equivalent to (19) for $`x_i = q_i \forall i`$.

# TBR

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
