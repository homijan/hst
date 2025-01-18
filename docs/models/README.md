# HST

## Standard multi-resolution analysis

Assuming $`\gamma_j = 1`$, the multi-resolution decomposition from *fine data* $`\phi_0`$ to *coarse data* $`\phi_J`$

$`\phi_j = G_{j-1} \phi_{j-1}\quad and\quad \overline{\phi}_j = \overline{G}_{j-1} \phi_{j-1},~\text{(4, Ref3)}`$

can be *inverted* to reconstruct the original *fine data* $`\phi_0`$ from multi-resolution **low & high frequency components**
$`[\phi_J, \overline{\phi}_J, .., \overline{\phi}_1]`$ as

$`G_{j-1}^H \phi_j = G_{j-1}^H G_{j-1} \phi_{j-1} \quad and\quad \overline{G}_{j-1}^H \overline{\phi}_j = \overline{G}_{j-1}^H \overline{G}_{j-1} \phi_{j-1}`$

$`\Rightarrow \quad \phi_{j-1} = G_{j-1}^H \phi_j + \overline{G}_{j-1}^H \overline{\phi}_j,~\text{(8, Ref3)}`$

which holds only if the **wavelet orthonormal filters** define a unitary transformation satisfying

$`\overline{G}_j G_j = G_j \overline{G}_j = 0\quad and\quad G_j^H G_j + \overline{G}_j^H \overline{G}_j = I,~\text{(7, Ref3)}`$

## Nonlinear multi-resolution analysis

Similar to low-frequency  multi-resolution decompostion in $`\text{(4, Ref3)}`$, define a nonlinear low-frequency decomposition

$`S_j = \rho\left(G_{j-1} S_{j-1}\right)\quad and\quad \overline{S}_j = \overline{\rho}\left(\overline{G}_{j-1} S_{j-1}\right),~(1)`$

where $\rho$ and $\overline{\rho}$ are invertible functions, $`\rho\left(\rho^{-1}(f)\right) = \rho^{-1}\left(\rho(f)\right) = f`$ and $`\overline{\rho}\left(\overline{\rho}^{-1}(f)\right) = \overline{\rho}^{-1}\left(\overline{\rho}(f)\right) = f`$, acting on the data vector element by element, and the recursion $`(1)`$ starts from the *fine data* $`\phi_0`$ as

$`S_0 = \rho(\phi_0).~(2)`$

Then, similar to $`\text{(8, Ref3)}`$, the decomposition $(1)$ can be *inverted* to reconstruct the original *fine data* $`\phi_0`$ from *nonlinear multi-resolution* **low & high frequency components**
$`[S_J, \overline{S}_J, .., \overline{S}_1]`$ as

$`G_{j-1}^H\rho^{-1}\left(S_j\right) = G_{j-1}^H G_{j-1} S_{j-1}\quad and\quad \overline{G}_{j-1}^H \overline{\rho}^{-1}\left(\overline{S}_j\right) = \overline{G}_{j-1}^H\overline{G}_{j-1} S_{j-1}`$

$`\Rightarrow\quad S_{j-1} = G_{j-1}^H\rho^{-1}\left(S_j\right) + \overline{G}_{j-1}^H \overline{\rho}^{-1}\left(\overline{S}_j\right),~(3)`$

where $(3)$ holds because $`\rho`$ and $`\overline{\rho}`$ are invertible and wavelet filters $G_j$ and $\overline{G}_j$ define a unitary transformation $`\text{(7, Ref3)}`$.

Note, that the recursive application of $`(3)`$ leads to a reconstruction of the *fine data* vector $`S_0`$. This is simply followed by the inversion of $`(2)`$

$`\phi_0 = g^{-1}\left( S_0 \right),~(4)`$

thus completing the *nonlinear multi-resolution decomposition and reconstruction*!

# References
1. `MST_Javier_Minguillon.pdf`
2. M. Glinsky, *A transformational approach to collective behavior*, `arXiv:2410.08558v4` (2025)
3. T. Marchand et al, *Wavelet Conditional Renormalization Group*, `arXiv:2207.04941v1` (2022)

# Notes

$`\left[ W_k f\right](x) = \left[\psi_k * f \right](x) = \int \psi_k(x^\prime) f(x - x^\prime) dx^\prime`$

$`R_0(z) = z + \exp^{i \arg (z)}`$

$`\rho(z) = i \ln (R_0(z))`$

$`S_m(x) = \rho\left( \left[ W_{k_m} S_{m-1} \right](x)\right) = \rho\left(\left[ \psi_{k_m} * S_{m-1} \right](x)\right)`$

$`S_0(x) = \rho(f(x))`$

$`S_1(x) = \rho\left(\left[ \psi_{k_1} * S_0 \right](x)\right) = \rho\left(\left[\psi_{k_1} * \rho(f)\right](x)\right)`$

$`S_2(x) = \rho\left(\left[ \psi_{k_2} * S_1 \right](x)\right) = \rho\left(\left[\psi_{k_2} * \rho\left(\left[\psi_{k_1} * \rho(f)\right]\right)\right](x)\right)`$

$`S_3(x) = \rho\left(\left[ \psi_{k_2} * S_2 \right](x)\right) = \rho\left(\left[\psi_{k_3} * \rho\left(\left[\psi_{k_2} * \rho\left(\left[\psi_{k_1} * \rho(f)\right]\right)\right]\right)\right](x)\right)`$
