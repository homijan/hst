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

$`\phi_j = G_{j-1} \phi_{j-1}\quad and\quad \overline{\phi}_j = \overline{G}_{j-1} \phi_{j-1}~\text{(4, Ref3)}`$

can be *inverted* to reconstruct the original *fine data* from multi-resolution **low & high frequency components**

$`G_{j-1}^T \phi_j = G_{j-1}^T G_{j-1} \phi_{j-1} \quad and\quad \overline{G}_{j-1}^T \overline{\phi}_j = \overline{G}_{j-1}^T \overline{G}_{j-1} \phi_{j-1}\quad \Rightarrow \quad \phi_{j-1 = }`$

# References
1. `MST_Javier_Minguillon.pdf`
2. M. Glinsky, *A transformational approach to collective behavior*, `arXiv:2410.08558v4` (2025)
3. T. Marchand et al, *Wavelet Conditional Renormalization Group*, `arXiv:2207.04941v1` (2022)
