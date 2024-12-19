# Heisenberg scattering transform (HST)

## Install

Our `hst` python package can be installed locally following these steps

1. Clone the package source code (now public) `git clone https://github.com/homijan/hst.git`

2. Install the python package locally (editable mode) `pip install -e hst`

3. Check if you can load the following from the module in `python` in terminal

   `>>> from hst.mrw1d import generate_G_operators, save_G_operators, data_decomposition, data_reconstruction`

## Demo

Jupyter notebook `hst/docs/demo_hst.ipynb` shows how to *decompose* (transform) and *reconstruct* (inverse transform) given randomly generated super-Gaussian pulses in 1D with `db1` or `db2` Daubechies wavelets. Set `verify_Gs = True` to see that `G_operators` correspodning to *mother* (`G` or `G_lo`) and *daughter* (`bar_G` or `G_hi`) wavelet operators are orthogonal

$`G \bar{G}^{\dagger} = \bar{G} G^{\dagger} = 0`$

and invertible

$`G^{\dagger}G + \bar{G}^{\dagger}\bar{G} = I.`$

The matrix `input_data` can be adjusted using any 1D dataset that has binate structure corresponding to multiresoltion levels `Nlevels`. 

## Multiresolution wavelet methodology

We follow the theory of multiresolution analysis described in "Marchand et al, *Wavelet Conditional Renormalization Group* (2022)"
visualized in Fig. 2, more precisely applied as

$`\phi_j = \gamma_j^{-1} G_{j-1} \phi_{j-1} \quad\text{and}\quad \bar{\phi}_j = \gamma_j^{-1} \bar{G}_{j-1} \phi_{j-1},~(4)`$

and

$`\phi_{j-1} = \gamma_j G_{j-1}^\dagger \phi_j + \gamma_j \bar{G}_{j-1}^\dagger \bar{\phi}_j.~(8)`$
## Daubechies Wavelets
- Discrete wavelets defined through coefficients $`a_k \neq 0`$ for handful of $`k`$s
- Orthogonality $`b_k = (-1)^k a^*_{1-k},~(1)`$ (holds for any discrete wavelet)
- The orthogonal low-frequency operator $G$ and high frequency $\bar{G}$ are contructed in terms of discrete convolution through $`a_k`$s and $`b_k`$s and by binating the resolution of the data  

### Real DAUB4 `db2`
- `dec_lo` filter coefficients $`[a_{-1}, a_0, a_1, a_2]`$ = `[-0.12940952+0.j, 0.22414387+0.j, 0.8365163 +0.j, 0.48296291+0.j]`
- corresponding `dec_hi` filter coefficients $`[b_{-1}, b_0, b_1, b_2] = [-a_{2}, a_1, -a_0, a_{-1}]`$ = `[-0.48296291+0.j  0.8365163 +0.j -0.22414387+0.j -0.12940952-0.j]`
- a simple **one-level** 4x8 $G$ and $\bar{G}$ operators with *outise-zero-BC* (*right-to-left* reversed to match `pywt`)

  $`G = \begin{bmatrix}
  \alpha a_1 & a_0 & a_{-1} & 0 & 0 & 0 & 0 & 0
  \\
  0 & a_2 & a_1 & a_0 & a_{-1} & 0 & 0 & 0
  \\
  0 & 0 & 0 & a_2 & a_1 & a_0 & a_{-1} & 0
  \\
  0 & 0 & 0 & 0 & 0 & a_2 & a_1 & \alpha a_0
  \end{bmatrix}`$

  $`\bar{G} = \begin{bmatrix}
  -\alpha a_0 & a_1 & -a_2 & 0 & 0 & 0 & 0 & 0
  \\
  0 & a_{-1} & -a_0 & a_1 & -a_2 & 0 & 0 & 0
  \\
  0 & 0 & 0 & a_{-1} & -a_0 & a_{1} & -a_2 & 0
  \\
  0 & 0 & 0 & 0 & 0 & a_{-1} & -a_0 & \alpha a_1
  \end{bmatrix}`$

  where $`\alpha`$ scales $a$ coeficients adjusting the BC to maintain $G$ and $\bar{G}$ orthogonal and invertible. 

- BC scaling $\alpha$ is obtained from $`G\bar{G}^H = 0`$ and $`G^{H}G + \bar{G}^{H} \bar{G} = I`$ leading to

  $`\begin{align}
  -\alpha a_1 (\alpha a_0)^* - a_0 a_1^* + a_{-1} a_2^* &= 0,
  \\
  (\alpha a_1)^* \alpha a_1 + (\alpha a_0)^* \alpha a_0 &= 1,
  \end{align}`$
  
  $`\begin{align}
  a_{-1} + a_0 + a_1 + a_2 = \sqrt{2},
  \\
  -a_{2} + a_1 -a_0 + a_{-1} = 0,
  \end{align}`$
  
### Complex Symmetric Daubechies Wavelets (SDW)

- Symmetry $`a_k = a_{1-k}`$
- In the case of **SDW2**, $`a_k \neq 0`$ for $`k = -2, -1, 0, 1, 2`$ and the SDW2 filter coeffs $`a_1 = 0.662912+0.171163j, a_2 = 0.110485-0.085581j, a_3 = -0.066291-0.085581j`$ [Lima, *Image Processing with Complex Daubechies Wavelets* (1997)],
- a simple **one-level** 4x8 $G$ and $\bar{G}$ operators with *outise-zero-BC* (*right-to-left* reversed to match `pywt`)

  $`G = \begin{bmatrix}
  \alpha a_1 & \beta a_0 & \delta a_{-1} & a_{-2} & 0 & 0 & 0 & 0
  \\
  a_3 & a_2 & a_1 & a_0 & a_{-1} & a_{-2} & 0 & 0
  \\
  0 & 0 & a_3 & a_2 & a_1 & a_0 & a_{-1} & a_{-2}
  \\
  0 & 0 & 0 & 0 & a_{3} & \delta a_2 & \beta a_1 & \alpha a_0
  \end{bmatrix},~(2)`$

  $`\bar{G} = \begin{bmatrix}
  -\alpha a_0 & \beta a_1 & -\delta a_2 & a_{3} & 0 & 0 & 0 & 0
  \\
  -a_{-2} & a_{-1} & -a_0 & a_1 & -a_2 & a_3 & 0 & 0
  \\
  0 & 0 & -a_{-2} & a_{-1} & -a_0 & a_{1} & -a_2 & a_3
  \\
  0 & 0 & 0 & 0 & -a_{-2} & \delta a_{-1} & -\beta a_0 & \alpha a_1
  \end{bmatrix},~(3)`$

  where $`\alpha, \beta, \delta`$ are scaling coeficients adjusting the BC to maintain $G$ and $\bar{G}$ orthogonal and invertible. 

- BC coefficients $\tilde{a}$ are obtained from $`G\bar{G}^H = 0`$ and $`G^{H}G + \bar{G}^{H} \bar{G} = I`$ leading to

  $`\begin{align}
  -\alpha a_1 (\alpha a_0)^* + \beta a_0 (\beta a_1)^* - \delta a_{-1} (\delta a_2)^* + a_{-2} a_3^* &= 0
  \Rightarrow \beta^2 a_1 a_1^* - \delta^2 a_2 a_2^* = \alpha^2 a_1 a_1^* - a_3 a_3^* \\
  -\alpha a_1 (a_{-2})^* + \beta a_0 a_{-1}^* - \delta a_{-1} a_0^* + a_{-2} a_1^* &= 0
  \Rightarrow \alpha a_1 a_3^* - \beta a_1 a_1^* + \delta a_2 a_1^* - a_3 a_1^* = 0
  \\    
  (\alpha a_1)^* \alpha a_1 + a_3^* a_3 + (\alpha a_0)^* \alpha a_0 + a_{-2}^* a_{-2} &= 1 \Rightarrow \alpha^2 = \frac{1 - 2 a_3^* a_3}{2 a_1^* a_1} 
  \end{align}`$

$`\begin{align}
\alpha &= \sqrt{\frac{1 - 2 a_3^2}{2 a_1^2}}~(4)
\\
\beta^2 a_1^2 - \delta^2 a_2^2 &= \frac{1}{2} - 2 a_3^2 \quad\Rightarrow\quad \frac{\left( \delta a_2 a_1^* + \sqrt{\frac{1 - 2 a_3^2}{2 a_1^2}} a_1 a_3^* - a_3 a_1^* \right)^2}{a_1^2} - \delta^2 a_2^2 - \frac{1}{2} + 2 a_3^2 = 0
\\
\delta a_2 a_1^* - \beta a_1^2 &= a_3 a_1^* - \sqrt{\frac{1 - 2 a_3^2}{2 a_1^2}} a_1 a_3^* \quad\Rightarrow\quad \beta = \frac{\delta a_2 a_1^* + \sqrt{\frac{1 - 2 a_3^2}{2 a_1^2}} a_1 a_3^* - a_3 a_1^*}{a_1^2}~(5)
\end{align}`$

$`(c_0 \delta + c_1)^2 + c_2 \delta^2 + c_3 = 0,~(6)`$

Mother wavelet (or *scaling function* using $`a_k`$)

`dec_lo [-0.04687482-0.06051491j, 0.07812469-0.06051491j, 0.46874957+0.12103052j, 0.46874957+0.12103052j, 0.07812469-0.06051491j, -0.04687482-0.06051491j]`

Daughter wavelet (or *wavelet* using $`b_k`$)

`dec_hi [-0.04687482+0.06051491j, -0.07812469-0.06051491j, 0.46874957-0.12103052j, -0.46874957+0.12103052j, 0.07812469+0.06051491j, 0.04687482-0.06051491j]`

Demonstration of orthogonality 

`dec_lo.dec_hi = (3.469446951953614e-18+0j)`

Demonstration of invertibility

`dec_lo^*.dec_lo + dec_hi^*.dec_hi = (0.9999974786819997+0j)`

