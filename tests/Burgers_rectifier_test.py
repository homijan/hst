import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

# Logarithmic operation on high-frequencies
eps = 1e-10
c_nln = 1e-2

def R_0(z : complex, outside=True) -> complex:
    z_bar = 2.0 * z / np.pi
    z_bar = np.array(z_bar).astype(complex)
    
    value_plus = (z_bar + np.sqrt( z_bar * z_bar - 1 )) / 1j
    value_minus = (z_bar - np.sqrt( z_bar * z_bar - 1 )) / 1j
    
    if outside:
        value = np.where(np.real(z) >= 0, value_plus, value_minus) # MOD: use value_plus for positive real part
    else:
        value = np.where(np.real(z) >= 0, value_minus, value_plus)
        
    return value

def R_0_inv(z : complex) -> complex:
    value = np.pi * (z * 1j + 1/(z * 1j) ) / 4.0
    return value

def R_inv(z):
    value = R_0_inv(np.exp(np.array(z).astype(complex) / 1j))
    return value

def R(z, outside=True):
    value = 1j * np.log(R_0(np.array(z).astype(complex), outside=outside))
    value = np.real(value) + 1j * np.maximum(np.imag(value), 0.0) # MOD: ensure non-negative imaginary part
    return value

def rho(z : complex):
    return z

def rho_inverse(z : complex):
    return z

def bar_rho(z : complex):
    return R(z)

def bar_rho_inverse(z : complex): 
    return R_inv(z)

# h inv test
if 1:
    def h_inv(z : complex):
        z_bar = 2.0 * z / np.pi
        z_bar = np.array(z_bar).astype(complex)
        
        value_plus = (z_bar + np.sqrt( z_bar * z_bar - 1 )) / 1j
        value_minus = (z_bar - np.sqrt( z_bar * z_bar - 1 )) / 1j
        if 0:
            value = np.where(np.abs(value_plus) > 1, value_plus, value_minus)
        else:
            value = np.where(np.real(z) >= 0, value_plus, value_minus)
        return value_plus, np.abs(value_plus), value_minus, np.abs(value_minus), value
    
    xx = np.linspace(-np.pi, np.pi, 200)
    hinv = h_inv(xx *1 - 0*np.pi/1000)
    plt.figure(figsize=(10, 5))
    plt.plot(xx, np.real(hinv[0]), label='real part +')
    plt.plot(xx, np.imag(hinv[0]), '*', label='imag part +')
    plt.plot(xx, np.abs(hinv[1]), label='abs +')
    plt.plot(xx, np.real(hinv[2]), 'o', label='real part -')
    plt.plot(xx, np.imag(hinv[2]), label='imag part -')
    plt.plot(xx, np.abs(hinv[3]), label='abs -')
    plt.plot(xx, np.real(hinv[4]), '--k', label='real part h inv')
    plt.plot(xx, np.imag(hinv[4]), ':k', label='imag part h inv')
    plt.plot(xx, np.abs(hinv[4]), '-.', label='abs h inv')
    plt.legend()
    plt.title('Inverse of h function')
    plt.xlabel('x')
    plt.grid()
    plt.show()
    
    exit()



# load Burgers data
input_data_file = 'data/Burgers_output.npz'
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

# plot the input data
if 0:
    plt.figure(figsize=(10, 5))
    plt.imshow(input_data, aspect='auto', origin='lower', cmap='viridis')
    plt.xlabel('t')
    plt.ylabel('x')
    plt.colorbar()
    plt.title('Input Data from Burgers Equation')
    #plt.show()

def test_function(x):
    """Test function to apply the rectifier"""
    return x
    return x**3/10

t_idx = -1
if 0:
    # Test function for rectifier
    def test_function(x):
        return (1j*np.abs(np.sin(x)) + np.cos(x))*1.5#*1.5 + 2 + 2j
    N = 150
    xx = np.linspace(-np.pi, np.pi, N)
    input_data = test_function(xx)
else:
    # select time index
    t_idx = 14
    xx = x
    input_data = input_data[:, t_idx] #+ 1j * input_data[:, t_idx][::-1]
    input_data = input_data #- 0.5#* 15 - 10


plt.figure(figsize=(22, 8))
n = 3
m = 5
for i in range(n):
    reMax = np.max(np.real(input_data))
    imMax = np.max(np.imag(input_data))
    maxMax = max(reMax, imMax)*1.2

    # Apply rectifier
    outside = True
    rect_data = R(input_data, outside=outside)
    rect0_data = R_0(input_data, outside=outside)
    #rect_data2 = R(input_data, outside=False)
    #rect_data = rect_data + rect_data2
    print(np.real(rect_data))
    print(np.imag(rect_data))

    reMax = np.max(np.real(rect_data))
    imMax = np.max(np.imag(rect_data))
    maxMax2 = max(reMax, imMax)*1.2

    reMax = np.max(np.real(rect0_data))
    imMax = np.max(np.imag(rect0_data))
    maxMax3 = max(reMax, imMax)*1.2
    
    # plot the input data
    plt.subplot(n, m, 1 + i*m)
    plt.plot(xx, np.real(input_data), label=f'real')
    plt.plot(xx, np.imag(input_data), label=f'imag')
    plt.title(f'Input data t={t_idx}')
    plt.xlabel('x')
    plt.legend()
    plt.subplot(n, m, 2 + i*m)
    plt.plot(xx, np.real(rect_data), label=f'real, outside={outside}')
    plt.plot(xx, np.imag(rect_data), label=f'imag, outside={outside}')
    #plt.plot(np.real(rect_data2), label=f'real, outside=F')
    #plt.plot(np.imag(rect_data2), label=f'imag, outside=F')
    plt.xlabel('x')
    plt.legend()
    plt.title(f'Rectified data R(z) t={t_idx}')
    plt.subplot(n, m, 3 + i*m)
    plt.scatter(np.real(input_data), np.imag(input_data), c=xx, cmap='viridis', s=1)
    plt.colorbar(label='x')
    plt.title('Input data Complex Plane')
    plt.xlim(-maxMax, maxMax)
    plt.ylim(-maxMax, maxMax)
    plt.xlabel('Re')
    plt.ylabel('Im')
    plt.subplot(n, m, 4 + i*m)
    plt.scatter(np.real(rect0_data), np.imag(rect0_data), c=xx, cmap='viridis', s=1)
    plt.colorbar(label='x')
    plt.title('R_0 Complex Plane')
    plt.xlim(-maxMax3, maxMax3)
    plt.ylim(-maxMax3, maxMax3)
    plt.xlabel('Re')
    plt.ylabel('Im')
    plt.subplot(n, m, 5 + i*m)
    plt.scatter(np.real(rect_data), np.imag(rect_data), c=xx, cmap='viridis', s=1)
    plt.xlim(-maxMax2, maxMax2)
    plt.ylim(-maxMax2, maxMax2)
    plt.xlabel('Re')
    plt.ylabel('Im')
    plt.colorbar(label='x')
    plt.title('Rectified data Complex Plane')
    
    input_data = rect_data  # use rectified data for the next iteration
    #input_data = np.real(rect_data)
    
plt.tight_layout()

if 0:
    spikes = np.zeros(len(rect_data))

    # compute positions of spikes
    for i in range(0, len(rect_data)):
        avg = (np.real(rect_data) + np.real(rect_data2)) / 2
        if not (np.real(rect_data[i]) < avg[i] and np.real(rect_data2[i]) > avg[i]): #or (ddx2 > 1 and ddx2 < 2)
            spikes[i] = 1
            plt.subplot(1, 2, 1)
            plt.plot(xx[i], test_function(xx[i]), 'ro')  # mark spike positions
            plt.subplot(1, 2, 2)
            plt.plot(i, 0, 'ro')  # mark spike positions

    print(f'Spike positions: {xx[spikes > 0]}')
    print(f'Spike f vals: {input_data[spikes > 0]}')
    print(f'Spikes count: {np.count_nonzero(spikes)}, fraction: {np.count_nonzero(spikes) / len(spikes)}')
    
    spike_idxs = np.where(spikes > 0)[0]
    
    idx = spike_idxs[0]
    print(f'Spike idx: {idx}')
    
    z_spike = np.complex64(2.0 * input_data[idx-1] / np.pi)
    value_plus = (z_spike + np.sqrt( z_spike * z_spike - 1 )) / 1j
    value_minus = (z_spike - np.sqrt( z_spike * z_spike - 1 )) / 1j    
    print(f'z_spike-1: {z_spike}, value_plus: {value_plus}, mag: {np.abs(value_plus):.30f}, +/-: {np.abs(value_plus) > 1.0}, value_minus: {value_minus}, {z_spike * z_spike - 1}')
    
    z_spike = np.complex64(2.0 * input_data[idx] / np.pi)
    value_plus = (z_spike + np.sqrt( z_spike * z_spike - 1 )) / 1j
    value_minus = (z_spike - np.sqrt( z_spike * z_spike - 1 )) / 1j
    print(f'z_spike: {z_spike}, value_plus: {value_plus}, mag: {np.abs(value_plus):.30f}, +/-: {np.abs(value_plus) > 1.0}, value_minus: {value_minus}, {z_spike * z_spike - 1}')
    
    z_spike = np.complex64(2.0 * input_data[idx+1] / np.pi)
    value_plus = (z_spike + np.sqrt( z_spike * z_spike - 1 )) / 1j
    value_minus = (z_spike - np.sqrt( z_spike * z_spike - 1 )) / 1j
    print(f'z_spike+1: {z_spike}, value_plus: {value_plus}, mag: {np.abs(value_plus):.30f}, +/-: {np.abs(value_plus) > 1.0}, value_minus: {value_minus}, {z_spike * z_spike - 1}')
    
    z_spike = np.complex64(2.0 * input_data[idx+2] / np.pi)
    value_plus = (z_spike + np.sqrt( z_spike * z_spike - 1 )) / 1j
    value_minus = (z_spike - np.sqrt( z_spike * z_spike - 1 )) / 1j
    print(f'z_spike+2: {z_spike}, value_plus: {value_plus}, mag: {np.abs(value_plus):.30f}, +/-: {np.abs(value_plus) > 1.0}, value_minus: {value_minus}, {z_spike * z_spike - 1}')
    
    
plt.show()

