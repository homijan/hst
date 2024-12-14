import numpy as np
import matplotlib.pyplot as plt
from hst.mrw1d import generate_G_operators, save_G_operators, data_decomposition, data_reconstruction

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
verify_Gs = False#True

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
#input_data = np.random.rand(data_length, Ndata)

print(f'Nlevels {Nlevels}, data_length {data_length}, Ndata {Ndata}')

# Data decomposition into wavelet basis
wavelets = ['db1', 'db2']
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
    # Generate data decomposition into (phi_J, bar_phi_J, .., bar_phi_1)
    print('Compute data decomposition:')
    decomposition = data_decomposition(G_operators, input_data, verify_Gs)
    print('Decomposition done!')
    # Coarse phi_J interpretation
    cA = decomposition[0]
    #print(f'wavelet {wavelet}, cA')
    #print(f'{cA}')

    # Generate data reconstruction from (phi_J, bar_phi_J, .., bar_phi_1) 
    print('Compute data reconstruction:')
    reconstructed_data = data_reconstruction(decomposition, G_operators)
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

# DAUB4 or db2
dec_lo = np.array([-0.12940952255126037, 0.2241438680420134, 0.8365163037378079, 0.48296291314453416])
# Symmetric (a_k = a_{1-k}) SWD2 with indexes [-2, -1, 0, 1, 2, 3]
a1 = 0.662912 + 0.171163j
a2 = 0.110485 - 0.085581j
a3 = -0.066291 - 0.085581j
# Symmetry and unitary scaling 1/sqrt(2)
dec_lo = np.array([a3 / 2.0**0.5, a2 / 2.0**0.5, a1 / 2.0**0.5, a1 / 2.0**0.5, a2 / 2.0**0.5, a3 / 2.0**0.5])

N = len(dec_lo)
dec_hi = np.zeros(N, dtype=complex)
for index in range(N):
    # offeset of the local k from python i index as k_loc = i + offset
    offset = int(1 - N / 2)
    k_loc = index + offset
    print(f'array index {index}, k_loc {k_loc}')
    # b_k = (-1)^k a^*_{1-k}
    dec_hi[index] = (-1)**k_loc * dec_lo[1 - k_loc - offset].conjugate()
    
print('SDW2 filters')
print(f'dec_lo {dec_lo}')
print(f'dec_hi {dec_hi}')
print(f'Orthogonality')
print(f'dec_lo.dec_hi {sum(dec_lo * dec_hi)}')
print('Invertibility')
print(f'dec_lo^*.dec_lo + dec_hi^*.dec_hi = {sum(dec_lo.conjugate() * dec_lo) + sum(dec_hi.conjugate() * dec_hi)}')

# Check orthogonality and invertibility with zero-BC
dec_lo_BC = dec_lo[2:]
dec_hi_BC = dec_hi[2:]
print(f'Orthogonality')
print(f'dec_lo_BC.dec_hi_BC {sum(dec_lo_BC * dec_hi_BC)}')
print('Invertibility')
print(f'dec_lo_BC^*.dec_lo_BC + dec_hi_BC^*.dec_hi_BC = {sum(dec_lo_BC.conjugate() * dec_lo_BC) + sum(dec_hi_BC.conjugate() * dec_hi_BC)}')


