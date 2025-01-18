# HST modeling

$`\left[ W_k f\right](x) = \left[\psi_k * f \right](x) = \int \psi_k(x^\prime) f(x - x^\prime) dx^\prime`$

$`R_0(z) = z + \exp^{i \arg (z)}`$

$`\rho(z) = i \ln (R_0(z))`$

$`S_m(x) = \rho\left( \left[ W_{k_m} S_{m-1} \right](x)\right) = \rho\left(\left[ \psi_{k_m} * S_{m-1} \right](x)\right)`$

$`S_0(x) = \rho(f(x))`$

$`S_1(x) = \rho\left(\left[ \psi_{k_1} * S_0 \right](x)\right) = \rho\left(\left[\psi_{k_1} * \rho(f)\right](x)\right)`$

$`S_2(x) = \rho\left(\left[ \psi_{k_2} * S_1 \right](x)\right) = \rho\left(\left[\psi_{k_2} * \rho\left(\left[\psi_{k_1} * \rho(f)\right]\right)\right](x)\right)`$

$`S_3(x) = \rho\left(\left[ \psi_{k_2} * S_2 \right](x)\right) = \rho\left(\left[\psi_{k_3} * \rho\left(\left[\psi_{k_2} * \rho\left(\left[\psi_{k_1} * \rho(f)\right]\right)\right]\right)\right](x)\right)`$

## Standard multi-resolution analysis

Assuming $`\gamma_j = 1`$, the multi-resolution decomposition from *fine data* $`\phi_0`$ to *coarse data* $`\phi_J`$

$`\phi_j = G_{j-1} \phi_{j-1}\quad and\quad \overline{\phi}_j = \overline{G}_{j-1} \phi_{j-1},~\text{(4, Ref3)}`$

can be *inverted* to reconstruct the original *fine data* from multi-resolution **low & high frequency components**
$`[\phi_J, \overline{\phi}_J, .., \overline{\phi}_1]`$ as

$`G_{j-1}^T \phi_j = G_{j-1}^T G_{j-1} \phi_{j-1} \quad and\quad \overline{G}_{j-1}^T \overline{\phi}_j = \overline{G}_{j-1}^T \overline{G}_{j-1} \phi_{j-1}\quad \Rightarrow \quad \phi_{j-1} = G_{j-1}^T \phi_j + \overline{G}_{j-1}^T \overline{\phi}_j,~\text{(8, Ref3)}`$

which holds only if the **wavelet orthonormal filters** define a unitary transformation satisfying

$`\overline{G}_j G_j = G_j \overline{G}_j = 0\quad and\quad G_j^T G_j + \overline{G}_j^T \overline{G}_j = I,~\text{(7, Ref3)}`$

## Nonlinear multi-resolution analysis

Similar to low-frequency  multi-resolution decompostion in $`\text{(4, Ref3)}`$, define a nonlinear low-frequency decomposition

$`S_j = \rho(G_{j-1} S_{j-1})\quad and\quad \overline{S}_j = \overline{\rho}(\overline{G}_{j-1} S_{j-1}),~(1)`$

where $\rho$ is an invertible function acting on the data vector element by element, and the recursion $`(1)`$ starts of the *fine data* as

$`S_0 = \rho(f).~(2)`$

# References
1. `MST_Javier_Minguillon.pdf`
2. M. Glinsky, *A transformational approach to collective behavior*, `arXiv:2410.08558v4` (2025)
3. T. Marchand et al, *Wavelet Conditional Renormalization Group*, `arXiv:2207.04941v1` (2022)
