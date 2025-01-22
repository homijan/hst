import numpy as np
from hst.wavelet_operators import verify_G_operators

def nonlinear_data_decomposition(G_operators, data, fnln, bar_fnln, verify_Gs=False):
    """Implementation of the nonlinear wavelet decomposition (S_J, bar_S_J, .., bar_S_1)"""
    if (verify_Gs):
        # Verify orthogonality and invertibility of G_operators at all levels
        verify_G_operators(G_operators)

    # Execute upward decompostion from fine to coarse
    # following the blue procedure in Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    decomposition = []
    # Starting (finest) level data
    phi_0 = data
    # Apply the low-frequency nonlinearity to obtain S_0
    S = fnln(phi_0)
    # The G_operators levels need to be reversed upward
    for G_lo, G_hi in reversed(G_operators):
        # Project high frequency data vector
        bar_S = G_hi.dot(S)
        # Apply the high-frequency nonlinearity
        bar_S = bar_fnln(bar_S)
        decomposition.append(bar_S)
        print(f'G_lo.shape {G_lo.shape}, S.shape {S.shape}, S count {S.shape[0]*S.shape[1]}, G_lo count_nonzero {np.count_nonzero(G_lo.toarray())}')
        # current level low frequency data vector
        S = G_lo.dot(S)
        # Apply the low-frequency nonlinearity
        S = fnln(S)
    # Add S_J (coarsest level S)
    decomposition.append(S)

    # Construct downward decomposition from coarse to fine
    # as vector (S_J, bar_S_J, bar_S_J-1,.., bar_S_1)
    # following Eq. 5 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    decomposition.reverse()

    return decomposition


def nonlinear_data_reconstruction(decomposition, G_operators, fnln_inv, bar_fnln_inv):
    """Implementation of data reconstruction from the wavelet coeffs (S_J, bar_S_J, .., bar_S_1).""" 
    # Reconstruct by a downward cascade starting with S_J
    # following the red procedure in Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    S = decomposition[0] 
    for i in range(len(G_operators)):
        bar_S = decomposition[i+1]
        # Apply the inverse low-frequency and high-frequency nonlinearities
        S = fnln_inv(S)
        bar_S = bar_fnln_inv(bar_S)
        # Get this level wavelet filters (orthonormal operators)
        G_lo, G_hi = G_operators[i]
        # Reconstruct S_{j-1} from low-ferquency S_j and high-frequency bar_S_j
        S = G_lo.conjugate(False).transpose(copy=False).dot(S) + G_hi.conjugate(False).transpose(copy=False).dot(bar_S)
        print(f'G_lo.H.shape {G_lo.conjugate(False).transpose(copy=False).shape}, bar_S.shape {bar_S.shape}, bar_S count {bar_S.shape[0]*bar_S.shape[1]}, G_lo count_nonzero {np.count_nonzero(G_lo.toarray())}')
    # For clarity we highlight that final S is on the lowest (finest) level
    # following Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    S_0 = S

    # Apply the inverse nonlinearity to reconstruct the data in its original linear form
    phi_0 = fnln_inv(S_0)

    return phi_0
