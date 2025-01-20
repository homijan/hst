import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from hst.wavelet_operators import generate_G_operators
from hst.nonlinmultres import nonlinear_data_decomposition, nonlinear_data_reconstruction

# Unit operation on low-frequencies
def nonlinear_function(f):
    return f
def nonlinear_function_inverse(g):
    return g
# Logarithmic operation on high-frequencies
c_nln = 1e-1
def bar_nonlinear_function(f):
    return 1j*np.log(f + c_nln)
def bar_nonlinear_function_inverse(g): 
    return np.exp(-1j*g) - c_nln

##############
# Input data #
##############
input_data_file = 'Burgers_output.npz'
data = np.load(input_data_file)
X = data['X']
U = data['U']
# Create the PCHIP interpolator
interp = PchipInterpolator(X, U)
# Preprocess the input data to fit a "close binatable grid"
x = np.linspace(min(X), max(X), 1024)
input_data = interp(x)
# Input data counts
data_length = input_data.shape[0]
n_data = input_data.shape[1]
print(f'Input data from file {input_data_file}, data_length {data_length}, n_data {n_data}')

#######################################
# HST decompostion and reconstruction #
#######################################
# Symmetric Debauchy wavelet sdw2
wavelet = 'sdw2'

# Number of scatters
n_levels = 5

print(f'Using wavelet {wavelet} on {n_levels} scattering levels within the Heisenberg scattering transform.')

# Generate orthogonal G_lo (aka G) and G_hi (aka bar_G) operators
G_operators = generate_G_operators(wavelet, n_levels, data_length)
print('G_operators generated.')
# Generate data decomposition into (S_J, bar_S_J, .., bar_S_1)
print('Compute data decomposition:')
decomposition = nonlinear_data_decomposition(G_operators, input_data, nonlinear_function, bar_nonlinear_function)
print('Decomposition done!')
# Coarse S_J interpretation
S_J = decomposition[0]
bar_S_1 = decomposition[n_levels]
print(f'S_J.shape {S_J.shape}, bar_S_1.shape {bar_S_1.shape}')

# Generate data reconstruction from (S_J, bar_S_J, .., bar_S_1) 
print('Compute data reconstruction:')
reconstructed_data = nonlinear_data_reconstruction(decomposition, G_operators, nonlinear_function_inverse, bar_nonlinear_function_inverse)
print('Reconstruction done!')

# Verification of direct and inverse multiresolution decomposition and reconstruction
print(f'Verfication of the nonlinear multiresolution decomposition and its inverse reconstruction:')
print(f'input_data.shape {input_data.shape}, reconstructed_data.shape {reconstructed_data.shape}')
print('input_data[:, :] - reconstructed_data[:, :]')
print(f'{input_data[:, :] - reconstructed_data[:, :]}')

# Plot results
n_plots = 5 
for i in range(n_data):
    if i % int(n_data / n_plots) == 0:
        plt.plot(x, input_data[:, i].real)
        plt.plot(x, reconstructed_data[:, i].real, '-.')
plt.show()
