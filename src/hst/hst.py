import numpy as np
from hst.wavelet_operators import verify_G_operators

# Logarithmic operation on high-frequencies
eps = 1e-10
c_nln = 1e-2

def R0(f : complex) -> complex:
    return f + np.exp(1j*np.angle(f))
def R0_inverse(f : complex) -> complex:
    return f - np.exp(1j*np.angle(f))
def rho(f):
    return f
def rho_inverse(f):
    return f
def bar_rho(f):
    return 1j*np.log(R0(f))
def bar_rho_inverse(f): 
    return R0_inverse(np.exp(-1j*f))


def hts_data_decomposition(G_operators, data, keep_bar_Sm=True, verify_Gs=False):
    """Implementation of the nonlinear wavelet decomposition (S_J, bar_S_J, .., bar_S_1)"""
    if (verify_Gs):
        # Verify orthogonality and invertibility of G_operators at all levels
        verify_G_operators(G_operators)

    # Execute upward decompostion from fine to coarse
    # following the blue procedure in Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    decomposition = []
    # Starting (finest) level data
    f0 = data
    # Apply the low-frequency nonlinearity to obtain S_0
    bar_Sj = 1j*np.log(R0(f0))
    # The G_operators levels need to be reversed upward
    for G_lo, G_hi in reversed(G_operators):
        # Project high frequency data vector
        bar_Sjp1 = G_hi.dot(bar_Sj)
        # Apply the high-frequency nonlinearity
        bar_Sjp1 = bar_rho(bar_Sjp1)
        
        print(f'G_lo.shape {G_lo.shape}, S.shape {bar_Sj.shape}, S count {bar_Sj.shape[0]*bar_Sj.shape[1]}, G_lo count_nonzero {np.count_nonzero(G_lo.toarray())}')
        # Project low frequency data vector
        Sj = G_lo.dot(bar_Sj)
        # Apply the low-frequency nonlinearity
        Sj = rho(Sj)
        
        decomposition.append(Sj)
        
        # Go to next level
        bar_Sj = bar_Sjp1
        
    # Add S_J (coarsest level S)
    if keep_bar_Sm:
        decomposition.append(bar_Sj)

    # Construct downward decomposition from coarse to fine
    # as vector (S_J, bar_S_J, bar_S_J-1,.., bar_S_1)
    # following Eq. 5 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    decomposition.reverse()

    return decomposition


def hst_data_reconstruction(decomposition, G_operators, keep_bar_Sm=True):
    """Implementation of data reconstruction from the wavelet coeffs (S_J, bar_S_J, .., bar_S_1).""" 
    # check sizes
    if (keep_bar_Sm and len(decomposition) != len(G_operators) + 1) or ((not keep_bar_Sm) and len(decomposition) != len(G_operators)):
        print("Wrong sizes of decomposition and G_operators array! Maybe change value of keep_bar_Sm.")
        return
    
    # Reconstruct by a downward cascade starting with S_J
    # following the red procedure in Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    if keep_bar_Sm:
        bar_Sj = decomposition[0]
    else:
        Sm = decomposition[0] 
        bar_Sj = G_operators[0].H.dot(rho_inverse(Sm))
              
    for i in range(int(keep_bar_Sm),len(G_operators)):
        Sj = decomposition[i]
        # Apply the inverse low-frequency and high-frequency nonlinearities
        Sj = rho_inverse(Sj)
        bar_Sj = bar_rho_inverse(bar_Sj)
        # Get this level wavelet filters (orthonormal operators)
        G_lo, G_hi = G_operators[i]
        # Reconstruct bar_S_{j-1} from low-ferquency S_j and high-frequency bar_S_j
        bar_Sj = G_lo.H.dot(Sj) + G_hi.H.dot(bar_Sj)
        print(f'G_lo.H.shape {G_lo.H.shape}, bar_S.shape {bar_Sj.shape}, bar_S count {bar_Sj.shape[0]*bar_Sj.shape[1]}, G_lo count_nonzero {np.count_nonzero(G_lo.toarray())}')
    # For clarity we highlight that final S is on the lowest (finest) level
    # following Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    bar_S0 = bar_Sj

    # Apply the inverse nonlinearity to reconstruct the data in its original linear form
    f0 = bar_rho_inverse(bar_S0)

    return f0
