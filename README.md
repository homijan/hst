# Heisenberg scattering transform (HST)

## Install

Our `hst` python package can be installed locally following these steps

1. Clone the package source code (now public) `git clone https://github.com/homijan/hst.git`

2. Install the python package locally (editable mode) `pip install -e hst`

3. Check if you can load the following from the module in `python` in terminal

   `>>> from hst.mrw1d import generate_G_operators, save_G_operators, data_decomposition, data_reconstruction`

## Demo

Copy the **demo** script from the repository `hst/tests/test_wavelet_decomposition.py` which shows how to *decompose* (transform) and *reconstruct* (inverse transform) any given 1D data with `db1` or `db2` Daubechies wavelets. Set `verify_Gs = True` to see that `G_operators` correspodning to *mother* (`G` or `G_lo`) and *daughter* (`bar_G` or `G_hi`) wavelet operators are orthogonal and invertible. 

The matrix `input_data` can be adjusted using your 1D dataset. 

## Multiresolution wavelet methodology

TBD
