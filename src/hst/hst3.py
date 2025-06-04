import numpy as np
from hst.wavelet_operators import verify_G_operators

# Logarithmic operation on high-frequencies
eps = 1e-10
c_nln = 1e-2

def R_0(z : complex, outside=True) -> complex:
    z_bar = 2.0 * z / np.pi
    z_bar = np.array(z_bar).astype(complex)
    
    value_plus = (z_bar + np.sqrt( z_bar * z_bar - 1 )) / 1j
    value_minus = (z_bar - np.sqrt( z_bar * z_bar - 1 )) / 1j
    
    if outside:
        value = np.where(np.real(z) >= 0, value_plus, value_minus)
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
    #value = np.real(z) + 1j * np.abs(np.imag(z)) # ensure non-negative imaginary part
    value = 1j * np.log(R_0(np.array(z).astype(complex), outside=outside))
    value = np.real(value) + 1j * np.maximum(np.imag(value), 0.0) # ensure non-negative imaginary part
    return value

def rho(z : complex):
    return z #R(z)

def rho_inverse(z : complex):
    return z #R_inv(z)

def bar_rho(z : complex):
    return R(z)

def bar_rho_inverse(z : complex): 
    return R_inv(z)


def hst3_data_decomposition(G_operators, data, verify_Gs=False):
    """Implementation of the nonlinear wavelet decomposition (S_J, bar_S_J, .., bar_S_1)"""
    if (verify_Gs):
        # Verify orthogonality and invertibility of G_operators at all levels
        verify_G_operators(G_operators)

    # Execute upward decompostion from fine to coarse
    # following the blue procedure in Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    # Starting (finest) level data
    f0 = data
    # Apply the low-frequency nonlinearity to obtain S_0
    bar_Sj = bar_rho(f0)
    decomposition = [bar_Sj]
    decompositionFull = [[bar_Sj]]
    decompositionFullWavelet = [[f0]]
    
    # The G_operators levels need to be reversed upward
    for G_lo, G_hi in reversed(G_operators):
        
        decompositionNext = []
        decompositionNextWavelet = []
        for bar_Sj in decomposition:
            # Project high frequency data vector
            bar_Sjp1 = G_hi.dot(bar_Sj)
            # Apply the high-frequency nonlinearity
            decompositionNextWavelet.append(bar_Sjp1)
            bar_Sjp1 = bar_rho(bar_Sjp1)
            
            #print(f'G_lo.shape {G_lo.shape}, S.shape {bar_Sj.shape}, S count {bar_Sj.shape[0]*bar_Sj.shape[1]}, G_lo count_nonzero {np.count_nonzero(G_lo.toarray())}')
            # Project low frequency data vector
            Sj = G_lo.dot(bar_Sj)
            # Apply the low-frequency nonlinearity
            decompositionNextWavelet.append(Sj)
            Sj = rho(Sj)
            
            decompositionNext.append(Sj)
            decompositionNext.append(bar_Sjp1)
            
        decomposition = decompositionNext
        decompositionFull.append(decomposition)
        decompositionFullWavelet.append(decompositionNextWavelet)
        
    # Add S_J (coarsest level S)

    # Construct downward decomposition from coarse to fine
    # as vector (S_J, bar_S_J, bar_S_J-1,.., bar_S_1)
    # following Eq. 5 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    #decomposition.reverse()

    return decomposition, decompositionFull, decompositionFullWavelet


def hst3_data_reconstruction(decomposition, G_operators):
    """Implementation of data reconstruction from the wavelet coeffs (S_J, bar_S_J, .., bar_S_1).""" 
    
    # Reconstruct by a downward cascade starting with S_J
    # following the red procedure in Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)              
    for i in range(len(G_operators)):
        # Get this level wavelet filters (orthonormal operators)
        G_lo, G_hi = G_operators[i]
        decompositionNext = []
        
        for j in range(0, len(decomposition), 2):
            print(j)
            Sj = decomposition[j]
            bar_Sj = decomposition[j + 1]
            # Apply the inverse low-frequency and high-frequency nonlinearities
            Sj = rho_inverse(Sj)
            bar_Sj = bar_rho_inverse(bar_Sj)
            # Reconstruct bar_S_{j-1} from low-ferquency S_j and high-frequency bar_S_j
            bar_Sj = G_lo.conjugate(False).transpose(copy=False).dot(Sj) + G_hi.conjugate(False).transpose(copy=False).dot(bar_Sj)
            
            decompositionNext.append(bar_Sj)
        
        decomposition = decompositionNext
            
    # For clarity we highlight that final S is on the lowest (finest) level
    # following Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    bar_S0 = decomposition[0]

    # Apply the inverse nonlinearity to reconstruct the data in its original linear form
    f0 = bar_rho_inverse(bar_S0)

    return f0
