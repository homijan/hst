import numpy as np
from hst.mrw1d import wavelet_decomposition, data_reconstruction

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
    decomposition, G_operators = wavelet_decomposition(wavelet, Nlevels, mat_data, verify_Gs)
    # Coarse phi_J interpretation
    cA = decomposition[0]
    print(f'wavelet {wavelet}, cA {cA}')

    data = data_reconstruction(decomposition, G_operators, Nlevels)

    print(f'mat_data.shape {mat_data.shape}, data.shape {data.shape}')
    print('mat_data[:, 0]')
    print(f'{mat_data[:, 0]}')
    print('data[:, 0]')
    print(f'{data[:, 0]}')
    
