import pywt
from pywt import wavedec
import numpy as np
from hst.mrw1d import generate_G_operators, data_decomposition

# MAIN
def format_array(arr):
    return "[%s]" % ", ".join(["%.14f" % x for x in arr])
def wavelet_info(wavelet):
    print(wavelet)
    print('decomposition_low')
    print(format_array(wavelet.dec_lo))
    print('decomposition_high')
    print(format_array(wavelet.dec_hi))
    print('reconstruction_low')
    print(format_array(wavelet.rec_lo))
    print('reconstruction_high')
    print(format_array(wavelet.rec_hi))
    return wavelet.dec_lo, wavelet.dec_hi, wavelet.rec_lo, wavelet.rec_hi 

# Print out essential properties of G_operators
verify_Gs = True

print('Demonstration that we apply db2 filter properly. Only one level of convolution.')
print('Note that the BCs are different. We use zero values of points of conv-kernel outside the dataset.')
Ndata = 1
Nlevels = 1
data_length = 2**(Nlevels+3)
print(f'Nlevels {Nlevels}, data_length {data_length}')
mat_data = np.random.rand(Ndata, data_length)

wavelet = 'db2'
tmp_data = mat_data[0][1:-1]
coeffs = wavedec(np.array([tmp_data]), wavelet, level=Nlevels)
dec_lo, dec_hi, rec_lo, rec_hi = wavelet_info(pywt.Wavelet(wavelet))

print('wavelet {wavelet}, pywt.wavedec coeffs: A, D1, D2, ...')
for coeff in coeffs:
    print(f'length {len(coeff[0])}, coeff {coeff}')
cA = coeffs[0]

data = mat_data.T
# Generate orthogonal G_lo (aka G) and G_hi (aka bar_G) operators
G_operators = generate_G_operators(wavelet, Nlevels, data.shape[0])
# Generate data decomposition into (phi_J, bar_phi_J, .., bar_phi_1)
# This test is strictly for pure multi-resolution without the nonlinearity
nonlinear_function = False
decomposition = data_decomposition(G_operators, data, nonlinear_function, verify_Gs)
# Coarse phi_J interpretation
cA_ours = decomposition[0]

print('Note that 1:-1 indexes of our implementation match pywt values')
print('This shows we use the filter/Gops correctly, with a different BC though (on purpose)!')
print(f'cA') 
print(f'{cA}')
print(f'cA_ours')
print(f'{cA_ours}')
print(f'cA - cA_ours.T {cA - cA_ours.T}')


print('Demonstration of working with db1 (our implementation is cross-checked against pywt) and db2 (we use only our implementation with zero-outside BC)')
Ndata = 20
Nlevels = 3
data_length = 2**(Nlevels+2)
print(f'Nlevels {Nlevels}, data_length {data_length}')
mat_data = np.random.rand(Ndata, data_length)

wavelet = 'db1'
print(f'TEST: wavelet {wavelet}')
dec_lo, dec_hi, rec_lo, rec_hi = wavelet_info(pywt.Wavelet(wavelet))
coeffs = wavedec(mat_data, wavelet, level=Nlevels)

#print('wavelet {wavelet}, pywt.wavedec coeffs: A, D1, D2, ...')
#for coeff in coeffs:
#    print(f'length {len(coeff[0])}, coeff {coeff}')
cA = coeffs[0]

data = mat_data.T
# Generate orthogonal G_lo (aka G) and G_hi (aka bar_G) operators
G_operators = generate_G_operators(wavelet, Nlevels, data.shape[0])
# Generate data decomposition into (phi_J, bar_phi_J, .., bar_phi_1)
# This test is strictly for pure multi-resolution without the nonlinearity
nonlinear_function = False
decomposition = data_decomposition(G_operators, data, nonlinear_function, verify_Gs)
# Coarse phi_J interpretation
cA_ours = decomposition[0]
#print(f'wavelet {wavelet}, cA_ours {cA_ours}')
print(f'{wavelet} verification: cA - cA_ours^T') 
print(f'{cA - cA_ours.T}')

wavelet = 'db2'
dec_lo, dec_hi, rec_lo, rec_hi = wavelet_info(pywt.Wavelet(wavelet))
print(f'TEST: wavelet {wavelet}')
data = mat_data.T
# Generate orthogonal G_lo (aka G) and G_hi (aka bar_G) operators
G_operators = generate_G_operators(wavelet, Nlevels, data.shape[0])
# Generate data decomposition into (phi_J, bar_phi_J, .., bar_phi_1)
# This test is strictly for pure multi-resolution without the nonlinearity
nonlinear_function = False
decomposition = data_decomposition(G_operators, data, nonlinear_function, verify_Gs)
# Coarse phi_J interpretation
cA_ours = decomposition[0]
print(f'wavelet {wavelet}, cA_ours {cA_ours}')
