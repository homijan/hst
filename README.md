# Heisenberg scattering transform (HST)

## Install

Our `hst` python package can be installed locally following these steps

1. Clone the package source code (now public) `git clone https://github.com/homijan/hst.git`

2. Install the python package locally (editable mode) `pip install -e hst`

3. Check if you can load the following from the module in `python` in terminal

   `>>> from hst.mrw1d import generate_G_operators, save_G_operators, data_decomposition, data_reconstruction`

## Demo

Jupyter notebook `hst/docs/demo_hst.ipynb` shows how to *decompose* (transform) and *reconstruct* (inverse transform) given randomly generated super-Gaussian pulses in 1D with `db1` or `db2` Daubechies wavelets. Set `verify_Gs = True` to see that `G_operators` correspodning to *mother* (`G` or `G_lo`) and *daughter* (`bar_G` or `G_hi`) wavelet operators are orthogonal

$`G \bar{G}^{\dagger} = \bar{G} G^{\dagger} = 0 \quad\Leftrightarrow\quad b_k = (-1)^k a^*_{1-k},`$

and invertible

$`G^{\dagger}G + \bar{G}^{\dagger}\bar{G} = I.`$

The matrix `input_data` can be adjusted using any 1D dataset that has binate structure corresponding to multiresoltion levels `Nlevels`. 

## Multiresolution wavelet methodology

We follow the theory of multiresolution analysis described in "Marchand et al, *Wavelet Conditional Renormalization Group* (2022)"
visualized in Fig. 2, more precisely applied as

$`\phi_j = \gamma_j^{-1} G_{j-1} \phi_{j-1} \quad\text{and}\quad \bar{\phi}_j = \gamma_j^{-1} \bar{G}_{j-1} \phi_{j-1},~(4)`$

and

$`\phi_{j-1} = \gamma_j G_{j-1}^T \phi_j + \gamma_j \bar{G}_{j-1}^T \bar{\phi}_j.~(8)`$




