import numpy as np
import matplotlib.pyplot as plt
from hst.wavelet_operators import generate_G_operators, save_G_operators
from hst.nonlinmultres import nonlinear_data_decomposition, nonlinear_data_reconstruction

# Multi-coefficient super-Gaussian data generation
def superGaussian(x, cs): 
    data_length = len(x)
    Ndata = cs.shape[0]
    input_data = np.zeros((data_length, Ndata))
    for i in range(Ndata):
        input_data[:, i] = cs[i, 0]*np.exp(-((x-cs[i, 1])/cs[i, 2])**(2*cs[i, 3]))
    return input_data

# Save G_operators
save_Gs = False#True
# Print out essential properties of G_operators
verify_Gs = True

Nlevels = 7
Ndata = 15 # Number of randomized super-Gaussian (small number to plot)
data_length = 2**Nlevels * 100 # constructed as a proper  power of 2

#Nlevels = 3
#Ndata = 10
#data_length = 1024 # constructed as a proper  power of 2

# Generate data
domain_size = 1.0
x = np.linspace(-domain_size/2.0, domain_size/2.0, data_length)
# Random super-Gaussian c0*exp(-((x-c1)/c2)^(2*c3))
input_coeffs = np.random.rand(Ndata, 4)
# Rescale coeff values, so c0 in (0, 1), c1 in (-domain_size/2, domain_size/2), c2 in (0, domain_size/2), and c3 in {1, .., 10}
for i in range(Ndata):
    # c0 amplitude scaling is already correct
    input_coeffs[i, 0] = input_coeffs[i, 0]
    # c1 shift scaling
    input_coeffs[i, 1] = 0.8 * domain_size * (input_coeffs[i, 1] - 0.5)
    # c2 width scaling
    input_coeffs[i, 2] = 0.125 * domain_size * input_coeffs[i, 2]
    # c3 super-Gaussian power scaling
    input_coeffs[i, 3] = int(10 * input_coeffs[i, 3])
input_data = superGaussian(x, input_coeffs)

print(f'Nlevels {Nlevels}, data_length {data_length}, Ndata {Ndata}')

# Proof of concept nonlinear function
c_nln = 1e-1
power_nln = 1.0 / 2.0
def bar_nonlinear_function(f):
    return (f + c_nln)**power_nln
def bar_nonlinear_function_inverse(g): 
    return g**(1.0 / power_nln) - c_nln
    # A dummy test to show the sensitivity of the inverse function
    #return g**(1.0 / (0.9999 * power_nln)) - c_nln
def nonlinear_function(f):
    return f
def nonlinear_function_inverse(g):
    return g

# Data decomposition into wavelet basis
wavelets = ['db1', 'db2', 'sdw2']
for wavelet in wavelets:
    print(f'TEST: wavelet {wavelet}')
    # Generate orthogonal G_lo (aka G) and G_hi (aka bar_G) operators
    G_operators = generate_G_operators(wavelet, Nlevels, data_length)
    print('G_operators generated.')
    if (save_Gs):
        # Save generated G_operators using dictionary
        file_name = f'G_operators_{wavelet}.npz'
        save_G_operators(G_operators, file_name)
        print(f'G_operators saved into {file_name}.')
    # Generate data decomposition into (S_J, bar_S_J, .., bar_S_1)
    print('Compute data decomposition:')
    decomposition = nonlinear_data_decomposition(G_operators, input_data, nonlinear_function, bar_nonlinear_function, verify_Gs)
    print('Decomposition done!')
    # Coarse S_J interpretation
    cA = decomposition[0]
    #print(f'wavelet {wavelet}, cA')
    #print(f'{cA}')

    # Generate data reconstruction from (S_J, bar_S_J, .., bar_S_1) 
    print('Compute data reconstruction:')
    reconstructed_data = nonlinear_data_reconstruction(decomposition, G_operators, nonlinear_function_inverse, bar_nonlinear_function_inverse)
    print('Reconstruction done!')

    # Verification of direct and inverse multiresolution decomposition and reconstruction
    print(f'Verfication of the multiresolution decomposition and inverse reconstruction:')
    print(f'input_data.shape {input_data.shape}, reconstructed_data.shape {reconstructed_data.shape}')
    print('input_data[:, :] - reconstructed_data[:, :]')
    print(f'{input_data[:, :] - reconstructed_data[:, :]}')

    # Plot results
    for i in range(Ndata):
        plt.plot(x, input_data[:, i])
        plt.plot(x, reconstructed_data[:, i], '-.')
    plt.show()
