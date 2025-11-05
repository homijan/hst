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
\pi(a, \tau) = \partial_{\dot{f}} \mathcal{L}_f^{if}(a, \tau) = \rho_0(a) \dot{f}(a, \tau).~(L6) 
\end{equation}`$

Finally, we obtain the ideal fluid action density $`\mathcal{S}^{if}_f(a, \tau)`$ for a given function $`f`$ from $`(L3)`$ using ideal fluid Lagrangian density $`(L5)`$ by solving an ordinary differential equation

$`\begin{equation}
\frac{d \mathcal{S}_f^{if}}{d \tau}(a, \tau) = \rho_0(a) \left( \frac{1}{2}\dot{f}(a, \tau)^2 - \varepsilon\left( \frac{\rho_0(a)}{\partial_a f(a, \tau)} \right) \right).~(L7)
\end{equation}`$

**Remark:**
*It should be noted, that $`f(a, \tau)`$ is the evolving spatial position of a fluid element $`a`$ in time, $`f(a, \tau) = x(a, \tau)`$ in 1D*.

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

where we assume initial condition $`x(a, 0) = a`$.

We conclude this ideal fluid example by rewritting action density from $`(L7)`$ into laboratory Eulerian coordinates by using $`(L12)`$ as

$`\begin{equation}
\frac{d \mathcal{S}_f^{if}}{d \tau}(x, \tau) = \rho_0(a(x, \tau)) \left( \frac{1}{2}v(x, \tau)^2 - \varepsilon(x, \tau) \right),~(L13)
\end{equation}`$

and conjugate momentum $`(L6)`$ in lab coordinates as

$`\begin{equation}
\pi(x, \tau) = \rho_0(a(x, t)) v(x, \tau).~(L14) 
\end{equation}`$

The action $`S^{if}[f](\tau)`$ can be obtained from the solution of $`(L7)`$ via $`(L4)`$.

Note, that, equivalently, can be obtained from solution of $`(L13)`$ via

$`\begin{equation}
S[f](\tilde{\tau}) = \int_\Omega \mathcal{S}_f^{if}(x, \tilde{\tau}) \frac{1}{J}d^3x,~(L15)
\end{equation}`$

where the integration is carried out over laboratory coordinates (not the fluid coordinates, requiring Jacobian scaling $`J = \partial_a x(a(x, \tau), \tau)`$).

# Discrete formulation

Let's have a mesh in laboratory coordinates at $`t=0`$ given by nodes $x_i$, then we define the Lagrangian *fluid* coordinate by $`a_i = x_i \forall i`$.

The discrete analog of the Jacobian reads
$`J_i = \frac{x_{i+1} - x_i}{a_{i+1} - a_i} = \frac{x_{i+1} - x_i}{\frac{m_{i+1} - m_i}{{(\rho_0)}_i}}`$

The correctness of the fluid coordinate to lab coordinates transformation can be checked via density and the following needs to hold

$`\begin{equation}
\rho(\frac{x_{i+1}(t) + x_{i}(t)}{2}) = \frac{\rho_0(\frac{a_{i+1}(t) + a_{i}(t)}{2})}{J_i}
\end{equation}`$
