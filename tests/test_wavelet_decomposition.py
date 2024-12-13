import numpy as np
from hst.mrw1d import generate_G_operators, data_decomposition, data_reconstruction

# Print out essential properties of G_operators
verify_Gs = False#True

Nlevels = 6
Ndata = 1000
data_length = 2**12 #1024 # should be a power of 2

# Generate data
mat_data = np.random.rand(data_length, Ndata)

print(f'Nlevels {Nlevels}, data_length {data_length}, Ndata {Ndata}')

# Data decomposition into wavelet basis
wavelets = ['db1', 'db2']
for wavelet in wavelets:
    print(f'TEST: wavelet {wavelet}')
    # Generate orthogonal G_lo (aka G) and G_hi (aka bar_G) operators
    G_operators = generate_G_operators(wavelet, Nlevels, data_length)
    print('G_operators generated.')
    # Save generated G_operators using dictionary
    file_name = f'G_operators_{wavelet}.npz'
    savez_dict = {}
    savez_dict['allow_pickle'] = True
    for index, level_Gops in enumerate(G_operators):
        key_lo = f'G_lo_{index}'
        key_hi = f'G_hi_{index}'
        print(f'saving key {key_lo}, level_Gops[0].shape {level_Gops[0].shape}, level_Gops[1].shape {level_Gops[1].shape}')
        print(f'saving key {key_hi}, level_Gops[1].shape {level_Gops[1].shape}')
        savez_dict[key_lo] = level_Gops[0]
        savez_dict[key_hi] = level_Gops[1]
    np.savez(file_name, **savez_dict)
    print(f'G_operators saved into {file_name}.')
    # Generate data decomposition into (phi_J, bar_phi_J, .., bar_phi_1)
    decomposition = data_decomposition(G_operators, mat_data, verify_Gs)
    print('Decomposition done!')
    # Coarse phi_J interpretation
    cA = decomposition[0]
    print(f'wavelet {wavelet}, cA {cA}')

    # Generate data reconstruction from (phi_J, bar_phi_J, .., bar_phi_1) 
    data = data_reconstruction(decomposition, G_operators)
    print('Reconstruction done!')

    # Verification of direct and inverse multiresolution decomposition and reconstruction
    print(f'Verfication of the multiresolution decomposition and inverse reconstruction:')
    print(f'mat_data.shape {mat_data.shape}, data.shape {data.shape}')
    print('mat_data[:, :] - data[:, :]')
    print(f'{mat_data[:, :] - data[:, :]}')
    
