import numpy as np
from hst.mrw1d import G_operators, data_decomposition, data_reconstruction

# Print out essential properties of G_operators
verify_Gs = True

Nlevels = 6
Ndata = 20
data_length = 1024 # Should be a power of 2

# Generate data
mat_data = np.random.rand(data_length, Ndata)

print(f'Nlevels {Nlevels}, data_length {data_length}, Ndata {Ndata}')

# Data decomposition into wavelet basis
wavelets = ['db1', 'db2']
for wavelet in wavelets:
    print(f'TEST: wavelet {wavelet}')
    # Generate orthogonal G_lo (aka G) and G_hi (aka bar_G) operators
    G_operators_upward = G_operators(wavelet, Nlevels, data_length)
    # Generate data decomposition into (phi_J, bar_phi_J, .., bar_phi_1)
    decomposition = data_decomposition(G_operators_upward, mat_data, verify_Gs)
    # Coarse phi_J interpretation
    cA = decomposition[0]
    print(f'wavelet {wavelet}, cA {cA}')

    # Generate data reconstruction from (phi_J, bar_phi_J, .., bar_phi_1) 
    data = data_reconstruction(decomposition, G_operators_upward)

    print(f'mat_data.shape {mat_data.shape}, data.shape {data.shape}')
    print('mat_data[:, 0]')
    print(f'{mat_data[:, 0]}')
    print('data[:, 0]')
    print(f'{data[:, 0]}')
    
