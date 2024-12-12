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
    cA = decomposition[len(decomposition)-1][0]
    print(f'wavelet {wavelet}, cA {cA}')

    # Construct vector (phi_J, bar_phi_J, bar_phi_J-1,.., bar_phi_1) as in
    phis = []
    phi_J = decomposition[len(decomposition)-1][0]
    phis.append(phi_J)
    for phi, bar_phi in reversed(decomposition):
        print(f'bar_phi.shape {bar_phi.shape}')
        phis.append(bar_phi)

    # Reverse the order of G_operators levels
    G_operators.reverse()

    data = data_reconstruction(phis, G_operators, Nlevels)
    
