import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from hst.wavelet_operators import generate_G_operators
from hst.nonlinmultres import nonlinear_data_decomposition, nonlinear_data_reconstruction

#######################################
# Michael Glinsky - conformal mapping #
#######################################
def R_0_inv(z):
    value = np.pi * (z * 1j + 1/(z * 1j) ) / 4.0
    return value

def R_0(z, outside=True):
    z_bar = 2.0 * z / np.pi
    z_bar = np.array(z_bar).astype(complex)
    
    value_plus = (z_bar + np.sqrt( z_bar * z_bar - 1 )) / 1j
    value_minus = (z_bar - np.sqrt( z_bar * z_bar - 1 )) / 1j
    
    if outside:
        value = np.where(np.abs(value_plus) > 1.0, value_plus, value_minus)
    else:
        value = np.where(np.abs(value_plus) > 1.0, value_minus, value_plus)
        
    return value

def R_inv(z):
    value = R_0_inv(np.exp(np.array(z).astype(complex) / 1j))
    return value

def R(z, outside=True):
    value = 1j * np.log( R_0(np.array(z).astype(complex), outside=outside))
    return value
######################################

# Unit operation on low-frequencies
def nonlinear_function(f):
    return f
def nonlinear_function_inverse(g):
    return g
def bar_nonlinear_function(f):
    return f
def bar_nonlinear_function_inverse(g):
    return g


#def h_MG(z):
#    return 0.5 * (z + 1.0 / z)
#def h_MG_inverse(z):
#    #return z + (z * z - 1)**0.5
#    #if z.real < 0 or abs(z) < 1.0:
#    #if abs(z) < 1.0 and z.real < 0.0:
#    if abs(z) < 1.0:
#        return z - (z * z - 1)**0.5
#    else:
#        return z + (z * z - 1)**0.5
#
#
#    #if (np.abs(z) >= 1.0):
#    #    return z + (z * z - 1)**0.5
#    #else:
#    #    return z - (z * z - 1)**0.5
#    #indexes = np.abs(z) <= 1
#    #h_inv = z + (z * z - 1)**0.5
#    #h_inv[indexes] = z[indexes] - (z[indexes] * z[indexes] - 1)**0.5
#    #return h_inv
#
#def bar_nonlinear_function(z):
#    return 2 * 1j / np.pi * np.log(h_MG_inverse(z) / 1j)
#def bar_nonlinear_function_inverse(z):
#    return h_MG(1j * np.exp(np.pi / 2.0 / 1j * z))
#
#for i in range(100):
#    re = 4.0 * (np.random.rand() - 0.5)
#    im = 4.0 * (np.random.rand() - 0.5)
#    z = re + 1j * im
#    #print(f'z {z:.3e}, abs(z) {np.abs(z):.3e}')
#
#    h_h_inv_z = h_MG(h_MG_inverse(z))
#    h_inv_h_z = h_MG_inverse(h_MG(z))
#    nln_nln_inv_z = bar_nonlinear_function(bar_nonlinear_function_inverse(z))
#    nln_inv_nln_z = bar_nonlinear_function_inverse(bar_nonlinear_function(z))
#    eps = 1e-10
#    if np.abs(z) < 1.0 and (np.abs(z - h_h_inv_z) > eps or np.abs(z - h_inv_h_z) > eps):  
#        print(f'i {i}')
#        print(f'z {z:.3e}')
#        print(f'abs(z) {np.abs(z):.3e}')
#        print(f'h_MG(h_MG_inverse(z)) {h_h_inv_z:.3e}')
#        print(f'h_MG_inverse(h_MG(z)) {h_inv_h_z:.3e}')
#    if False:#np.abs(z - nln_nln_inv_z) > eps or np.abs(z - nln_inv_nln_z) > eps:
#        print(f'bar_nonlinear_function(bar_nonlinear_function_inverse(z)) {nln_nln_inv_z:.3e}')
#        print(f'bar_nonlinear_function_inverse(bar_nonlinear_function(z)) {nln_inv_nln_z:.3e}')


# Logarithmic operation on high-frequencies
eps = 1e-10
c_nln = 1e-2
# Note that that the one shift is turned off
def R0(f):
    return f + c_nln #+ f / (abs(f) + eps)
def R0_inverse(g):
    return g - c_nln #- g / (abs(g) + eps)
def MG_bar_nonlinear_function(f):
    return 1j*np.log(R0(f))
def MG_bar_nonlinear_function_inverse(g): 
    return R0_inverse(np.exp(-1j*g))

# log scatter from S_0 to coarser structures
#def bar_nonlinear_function(f):
#    return MG_bar_nonlinear_function(f)
#def bar_nonlinear_function_inverse(g):
#    return MG_bar_nonlinear_function_inverse(g)

# log scatter from S_J to finer structures
#def bar_nonlinear_function(f):
#    return MG_bar_nonlinear_function_inverse(f)
#def bar_nonlinear_function_inverse(g):
#    return MG_bar_nonlinear_function(g) 


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
print('Compute data decomposition (S_J, bar_S_J, .., bar_S_1):')
decomposition = nonlinear_data_decomposition(G_operators, input_data, nonlinear_function, bar_nonlinear_function)
print('Decomposition done!')
# Coarse S_J interpretation
S_J = decomposition[0]
bar_S_J = decomposition[1]
bar_S_1 = decomposition[n_levels]
print(f'S_J.shape {S_J.shape}, bar_S_J.shape {bar_S_J.shape}, bar_S_1.shape {bar_S_1.shape}')

# Generate data reconstruction from (S_J, bar_S_J, .., bar_S_1) 
print('Compute data reconstruction from (S_J, bar_S_J, .., bar_S_1):')
reconstructed_data = nonlinear_data_reconstruction(decomposition, G_operators, nonlinear_function_inverse, bar_nonlinear_function_inverse)
print('Reconstruction done!')

# Verification of direct and inverse multiresolution decomposition and reconstruction
print(f'Verfication of the nonlinear multiresolution decomposition and its inverse reconstruction:')
print(f'input_data.shape {input_data.shape}, reconstructed_data.shape {reconstructed_data.shape}')
print('input_data[:, :] - reconstructed_data[:, :]')
print(f'{input_data[:, :] - reconstructed_data[:, :]}')

#################
# Visualization #
#################
# Prepare x vectors for each level by "bination"
x_levels = []
xlevel = x
for i in range(n_levels):
    xlevel = xlevel[::2]
    print(f'xlevel.shape {xlevel.shape}')
    # Add x coordinates for bar_S_i
    x_levels.append(xlevel)
# Reverse the order to make bar_S_J x coords first
x_levels.reverse()

# Plot results
n_plots = 5
dsnp = int(n_data / n_plots)
print(f'dsnp {dsnp}')
fig1, ax = plt.subplots() 
fig2, axs = plt.subplots(n_plots, 3)

i_plot = 0
for i in range(n_data):
    if i % dsnp == 0:
        print(f'Snapshot {i}')
        ax.plot(x, input_data[:, i].real)
        ax.plot(x, reconstructed_data[:, i].real, '-.')
        # S_J
        ax_R = axs[i_plot, 0].twinx()
        axs[i_plot, 0].plot(x_levels[0], decomposition[0][:, i].real, 'x-', label='real')
        ax_R.plot(x_levels[0], decomposition[0][:, i].imag, 'o-', label='imag')
        # bar_S_j
        for j in range(n_levels):
            # Skip the S_J
            component = j + 1
            axs[i_plot, 1].plot(x_levels[j], decomposition[component][:, i].real, 'x-', label=f'j={n_levels-j}')
            axs[i_plot, 2].plot(x_levels[j], decomposition[component][:, i].imag, 'o-', label=f'j={n_levels-j}') 
        i_plot = i_plot + 1
axs[0, 0].set_title(f'S_J={n_levels}')
axs[0, 0].legend()
axs[0, 1].set_title('Re(bar_S)')
axs[0, 1].legend()
axs[0, 2].set_title('Im(bar_S)')
axs[0, 2].legend()
plt.show()
