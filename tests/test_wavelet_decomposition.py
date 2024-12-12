import numpy as np
from hst.mrw1d import wavelet_decomposition

# Print out essential properties of G_operators
verify_Gs = True

Nlevels = 10
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
    cA = decomposition[len(decomposition)-1][0]
    print(f'wavelet {wavelet}, cA {cA}')
