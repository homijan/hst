# Heisenberg scattering transform (HST)

## Install

Our `hst` python package can be installed locally following these steps

1. Clone the package source code (now public) `git clone https://github.com/homijan/hst.git`

2. Install the python package locally (editable mode) `pip install -e hst`

3. Check if you can load the following from the module in `python` in terminal

   `>>> from hst.mrw1d import generate_G_operators, save_G_operators, data_decomposition, data_reconstruction`

## Demo

Jupyter notebook `hst/docs/demo_hst.ipynb` shows how to *decompose* (transform) and *reconstruct* (inverse transform) given randomly generated super-Gaussian pulses in 1D with `db1` or `db2` Daubechies wavelets. Set `verify_Gs = True` to see that `G_operators` correspodning to *mother* (`G` or `G_lo`) and *daughter* (`bar_G` or `G_hi`) wavelet operators are orthogonal

$`G \overline{G}^{\dagger} = \overline{G} G^{\dagger} = 0`$

and invertible

$`G^{\dagger}G + \overline{G}^{\dagger}\overline{G} = I.`$

The matrix `input_data` can be adjusted using any 1D dataset that has binate structure corresponding to multiresoltion levels `Nlevels`. 

## Multiresolution wavelet methodology

We follow the theory of multiresolution analysis described in "Marchand et al, *Wavelet Conditional Renormalization Group* (2022)"
visualized in Fig. 2, more precisely applied as

$`\varphi_j = \gamma_j^{-1} G_{j-1} \varphi_{j-1} \quad\text{and}\quad \overline{\varphi}_j = \gamma_j^{-1} \overline{G}_{j-1} \varphi_{j-1},~(4)`$

and

$`\varphi_{j-1} = \gamma_j G_{j-1}^\dagger \varphi_j + \gamma_j \overline{G}_{j-1}^\dagger \overline{\varphi}_j.~(8)`$
## Daubechies Wavelets
- Discrete wavelets defined through coefficients $`a_k \neq 0`$ for handful of $`k`$s
- Orthogonality that holds for any discrete wavelet $`b_k = (-1)^k a^*_{1-k},~(0)`$
- The orthogonal low-frequency operator $G$ and high frequency $\overline{G}$ are contructed in terms of discrete convolution through $`a_k`$s and $`b_k`$s and by binating the resolution of the data  

### Real DAUB4 `db2`
- `dec_lo` filter coefficients $`[a_{-1}, a_0, a_1, a_2]`$ = `[-0.12940952+0.j, 0.22414387+0.j, 0.8365163 +0.j, 0.48296291+0.j]`
- corresponding `dec_hi` filter coefficients $`[b_{-1}, b_0, b_1, b_2] = [-a_{2}^*, a_1^*, -a_0^*, a_{-1}^*]`$ = `[-0.48296291+0.j  0.8365163 +0.j -0.22414387+0.j -0.12940952-0.j]`
- a simple **one-level** 4x8 $G$ and $\overline{G}$ operators with *outise-zero-BC* (*right-to-left* reversed to match `pywt`)

  $`G = \begin{bmatrix}
  \alpha a_1 & a_0 & a_{-1} & 0 & 0 & 0 & 0 & 0
  \\
  0 & a_2 & a_1 & a_0 & a_{-1} & 0 & 0 & 0
  \\
  0 & 0 & 0 & a_2 & a_1 & a_0 & a_{-1} & 0
  \\
  0 & 0 & 0 & 0 & 0 & a_2 & a_1 & \alpha a_0
  \end{bmatrix}`$

  $`\overline{G} = \begin{bmatrix}
  -(\alpha a_0)^* & a_1^* & -a_2^* & 0 & 0 & 0 & 0 & 0
  \\
  0 & a_{-1}^* & -a_0^* & a_1^* & -a_2^* & 0 & 0 & 0
  \\
  0 & 0 & 0 & a_{-1}^* & -a_0^* & a_{1}^* & -a_2^* & 0
  \\
  0 & 0 & 0 & 0 & 0 & a_{-1}^* & -a_0^* & (\alpha a_1)^*
  \end{bmatrix}`$

  where $`\alpha`$ scales $a$ coeficients adjusting the BC to maintain $G$ and $\overline{G}$ orthogonal and invertible. 

- BC scaling $\alpha$ is obtained from $`G\overline{G}^H = 0`$ and $`G^{H}G + \overline{G}^{H} \overline{G} = I`$ leading to

  $`-\alpha a_1 (\alpha a_0)^* - a_0 a_1^* + a_{-1} a_2^* = 0,`$
  
### Complex Symmetric Daubechies Wavelets (SDW)

- Symmetry $`a_k = a_{1-k}`$
- In the case of **SDW2**, $`a_k \neq 0`$ for $`k = -2, -1, 0, 1, 2`$ and the SDW2 filter coeffs [Lima, *Image Processing with Complex Daubechies Wavelets* (1997)]

  $`a_1 = 0.662912+0.171163j, a_2 = 0.110485-0.085581j, a_3 = -0.066291-0.085581j,~(1)`$

  Mother wavelet (or *scaling function* using $`a_k`$)

  `dec_lo [-0.04687482-0.06051491j, 0.07812469-0.06051491j, 0.46874957+0.12103052j, 0.46874957+0.12103052j, 0.07812469-0.06051491j, -0.04687482-0.06051491j]`

  Daughter wavelet (or *wavelet* using $`b_k`$)

  `dec_hi [-0.04687482+0.06051491j, -0.07812469-0.06051491j, 0.46874957-0.12103052j, -0.46874957+0.12103052j, 0.07812469+0.06051491j, 0.04687482-0.06051491j]`

  Demonstration of orthogonality 

  `dec_lo.dec_hi = (3.469446951953614e-18+0j)`

  Demonstration of invertibility

  `dec_lo^*.dec_lo + dec_hi^*.dec_hi = (0.9999974786819997+0j)`

- a simple **one-level** 4x8 $G$ and $\overline{G}$ operators with *outise-zero-BC* (*right-to-left* reversed to match `pywt`) *succesful method!*

  $`G = \begin{bmatrix}
  \alpha a_1 & \beta a_0 & a_{-1} & a_{-2} & 0 & 0 & 0 & 0
  \\
  a_3 & a_2 & a_1 & a_0 & a_{-1} & a_{-2} & 0 & 0
  \\
  0 & 0 & a_3 & a_2 & a_1 & a_0 & a_{-1} & a_{-2}
  \\
  0 & 0 & 0 & 0 & a_{3} & a_2 & \beta a_1 & \alpha a_0
  \end{bmatrix} = \begin{bmatrix}
  \alpha a_1 & \beta a_1 & a_2 & a_3 & 0 & 0 & 0 & 0
  \\
  a_3 & a_2 & a_1 & a_1 & a_2 & a_3 & 0 & 0
  \\
  0 & 0 & a_3 & a_2 & a_1 & a_1 & a_2 & a_3
  \\
  0 & 0 & 0 & 0 & a_{3} & a_2 & \beta a_1 & \alpha a_1
  \end{bmatrix},~(2)`$

  $`\overline{G} = \begin{bmatrix}
  -(\alpha a_0)^* & (\beta a_1)^* & -a_2^* & a_3^* & 0 & 0 & 0 & 0
  \\
  -a_{-2}^* & a_{-1}^* & -a_0^* & a_1^* & -a_2^* & a_3^* & 0 & 0
  \\
  0 & 0 & -a_{-2}^* & a_{-1}^* & -a_0^* & a_{1}^* & -a_2^* & a_3^*
  \\
  0 & 0 & 0 & 0 & -a_{-2}^* & a_{-1}^* & -(\beta a_0)^* & (\alpha a_1)^*
  \end{bmatrix} = \begin{bmatrix}
  -(\alpha a_1)^* & (\beta a_1)^* & -a_2^* & a_3^* & 0 & 0 & 0 & 0
  \\
  -a_3^* & a_2^* & -a_1^* & a_1^* & -a_2^* & a_3^* & 0 & 0
  \\
  0 & 0 & -a_3^* & a_2^* & -a_1^* & a_1^* & -a_2^* & a_3^*
  \\
  0 & 0 & 0 & 0 & -a_3^* & a_2^* & -(\beta a_1)^* & (\alpha a_1)^*
  \end{bmatrix},~(3),`$

  **Solution to the orthogonality constraint on BC**

  `alpha = (a_2*(a_2/(a_2 + a_3) + a_3*sqrt((a_1 - a_2 - a_3)*(a_1 + a_2 + a_3))/(a_1*(a_2 + a_3))) - a_2 + a_3)/a_3 = (1.0204969605748353+0.015876976569495486j)`
  
  `beta = a_2/(a_2 + a_3) + a_3*sqrt((a_1 - a_2 - a_3)*(a_1 + a_2 + a_3))/(a_1*(a_2 + a_3)) = (1.012298185699968-0.015876906424713344j)`

Explicit form of the orthogonal operators $`G`$ and $`\overline{G}`$ for SDW2

`Minimalistic mat_G with BCs
[[ 0.67378213+0.18519636j  0.67378215+0.162743j    0.110485  -0.085581j
  -0.066291  -0.085581j    0.        +0.j          0.        +0.j
   0.        +0.j          0.        +0.j        ]
 [-0.066291  -0.085581j    0.110485  -0.085581j    0.662912  +0.171163j
   0.662912  +0.171163j    0.110485  -0.085581j   -0.066291  -0.085581j
   0.        +0.j          0.        +0.j        ]
 [ 0.        +0.j          0.        +0.j         -0.066291  -0.085581j
   0.110485  -0.085581j    0.662912  +0.171163j    0.662912  +0.171163j
   0.110485  -0.085581j   -0.066291  -0.085581j  ]
 [ 0.        +0.j          0.        +0.j          0.        +0.j
   0.        +0.j         -0.066291  -0.085581j    0.110485  -0.085581j
   0.67378215+0.162743j    0.67378213+0.18519636j]]
Minimalistic bar_mat_G with BCs
[[-0.67378213+0.18519636j  0.67378215-0.162743j   -0.110485  -0.085581j
  -0.066291  +0.085581j    0.        +0.j          0.        +0.j
   0.        +0.j          0.        +0.j        ]
 [ 0.066291  -0.085581j    0.110485  +0.085581j   -0.662912  +0.171163j
   0.662912  -0.171163j   -0.110485  -0.085581j   -0.066291  +0.085581j
   0.        +0.j          0.        +0.j        ]
 [ 0.        +0.j          0.        +0.j          0.066291  -0.085581j
   0.110485  +0.085581j   -0.662912  +0.171163j    0.662912  -0.171163j
  -0.110485  -0.085581j   -0.066291  +0.085581j  ]
 [ 0.        +0.j          0.        +0.j          0.        +0.j
   0.        +0.j          0.066291  -0.085581j    0.110485  +0.085581j
  -0.67378215+0.162743j    0.67378213-0.18519636j]]
Orthogonality mat_G.bar_mat_G^H
[[ 3.18321758e-16-1.49186219e-16j  3.46944695e-18-2.08166817e-17j
   0.00000000e+00+0.00000000e+00j  0.00000000e+00+0.00000000e+00j]
 [ 3.46944695e-18-2.08166817e-17j -1.04083409e-17+1.73472348e-17j
  -6.93889390e-18+0.00000000e+00j  0.00000000e+00+0.00000000e+00j]
 [ 0.00000000e+00+0.00000000e+00j -6.93889390e-18+0.00000000e+00j
  -1.04083409e-17+1.73472348e-17j -1.04083409e-17+1.73472348e-17j]
 [ 0.00000000e+00+0.00000000e+00j  0.00000000e+00+0.00000000e+00j
  -1.04083409e-17+1.73472348e-17j -3.29597460e-16+1.38777878e-16j]]
Invertibility mat_G^H.mat_G + bar_mat_G^H.bar_mat_G
[[ 9.99997310e-01+0.00000000e+00j  0.00000000e+00-1.64401975e-08j
   2.57823578e-07+0.00000000e+00j  0.00000000e+00+1.76532138e-08j
  -1.07148000e-07+0.00000000e+00j  0.00000000e+00+0.00000000e+00j
   0.00000000e+00+0.00000000e+00j  0.00000000e+00+0.00000000e+00j]
 [ 0.00000000e+00+1.64401975e-08j  9.99997437e-01+0.00000000e+00j
   0.00000000e+00+9.56577859e-09j  2.59911753e-07+0.00000000e+00j
   0.00000000e+00+0.00000000e+00j -1.07148000e-07+0.00000000e+00j
   0.00000000e+00+0.00000000e+00j  0.00000000e+00+0.00000000e+00j]
 [ 2.57823578e-07+0.00000000e+00j  0.00000000e+00-9.56577859e-09j
   9.99997479e-01+0.00000000e+00j  0.00000000e+00+0.00000000e+00j
   2.63044000e-07+0.00000000e+00j  0.00000000e+00+0.00000000e+00j
  -1.07148000e-07+0.00000000e+00j  0.00000000e+00+0.00000000e+00j]
 [ 0.00000000e+00-1.76532138e-08j  2.59911753e-07+0.00000000e+00j
   0.00000000e+00+0.00000000e+00j  9.99997479e-01+0.00000000e+00j
   0.00000000e+00+0.00000000e+00j  2.63044000e-07+0.00000000e+00j
   0.00000000e+00+0.00000000e+00j -1.07148000e-07+0.00000000e+00j]
 [-1.07148000e-07+0.00000000e+00j  0.00000000e+00+0.00000000e+00j
   2.63044000e-07+0.00000000e+00j  0.00000000e+00+0.00000000e+00j
   9.99997479e-01+0.00000000e+00j  0.00000000e+00+0.00000000e+00j
   2.59911753e-07+0.00000000e+00j  0.00000000e+00-1.76532138e-08j]
 [ 0.00000000e+00+0.00000000e+00j -1.07148000e-07+0.00000000e+00j
   0.00000000e+00+0.00000000e+00j  2.63044000e-07+0.00000000e+00j
   0.00000000e+00+0.00000000e+00j  9.99997479e-01+0.00000000e+00j
   0.00000000e+00-9.56577859e-09j  2.57823578e-07+0.00000000e+00j]
 [ 0.00000000e+00+0.00000000e+00j  0.00000000e+00+0.00000000e+00j
  -1.07148000e-07+0.00000000e+00j  0.00000000e+00+0.00000000e+00j
   2.59911753e-07+0.00000000e+00j  0.00000000e+00+9.56577859e-09j
   9.99997437e-01+0.00000000e+00j  0.00000000e+00+1.64401975e-08j]
 [ 0.00000000e+00+0.00000000e+00j  0.00000000e+00+0.00000000e+00j
   0.00000000e+00+0.00000000e+00j -1.07148000e-07+0.00000000e+00j
   0.00000000e+00+1.76532138e-08j  2.57823578e-07+0.00000000e+00j
   0.00000000e+00-1.64401975e-08j  9.99997310e-01+0.00000000e+00j]]`
