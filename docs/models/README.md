# HST modeling

$`\left[ W_k f\right](x) = \left[\psi_k * f \right](x) = \int \psi_k(x^\prime) f(x - x^\prime) dx^\prime`$

$`R_0(z) = z + \exp^{i \arg (z)}`$

$`\rho(z) = i \ln (R_0(z))`$

$`S_0(x) = \rho(f(x))`$

$`S_m(x) = \rho\left( \left[ W_{k_m} S_{m-1} \right](x)\right) = \rho\left(\left[ \psi_{k_m} * S_{m-1} \right](x)\right)`$

$`S_1(x) = \rho\left([\psi_{k_1} * \rho(f)](x)\right)`$

$`S_2(x) = \rho\left([\psi_{k_2} * \rho\left([\psi_{k_1} * \rho(f)]\right)](x)\right)`$

$`S_3(x) = \rho\left([\psi_{k_3} * \rho\left([\psi_{k_2} * \rho\left([\psi_{k_1} * \rho(f)]\right)]\right)](x)\right)`$

# References
1. `MST_Javier_Minguillon.pdf`
2. M. Glinsky, *A transformational approach to collective behavior*, `arXiv:2410.08558v4` (2025)
3. T. Marchand et al, *Wavelet Conditional Renormalization Group*, `arXiv:2207.04941v1` (2022)
