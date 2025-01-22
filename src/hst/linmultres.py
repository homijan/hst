import numpy as np
from hst.wavelet_operators import verify_G_operators

def linear_data_decomposition(G_operators, data, verify_Gs=False):
    """Implementation of the nonlinear wavelet decomposition (phi_J, bar_phi_J, .., bar_phi_1)"""
    if (verify_Gs):
        # Verify orthogonality and invertibility of G_operators at all levels
        verify_G_operators(G_operators)

    # Execute upward decompostion from fine to coarse
    # following the blue procedure in Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    decomposition = []
    # Starting (finest) level data
    phi = data
    # The G_operators levels need to be reversed upward
    for G_lo, G_hi in reversed(G_operators):
        # Project high frequency data vector
        bar_phi = G_hi.dot(phi)
        decomposition.append(bar_phi)
        print(f'G_lo.shape {G_lo.shape}, phi.shape {phi.shape}, phi count {phi.shape[0]*phi.shape[1]}, G_lo count_nonzero {np.count_nonzero(G_lo.toarray())}')
        # current level low frequency data vector
        phi = G_lo.dot(phi)
    # Add phi_J (coarsest level phi)
    decomposition.append(phi)

    # Construct downward decomposition from coarse to fine
    # as vector (phi_J, bar_phi_J, bar_phi_J-1,.., bar_phi_1)
    # following Eq. 5 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    decomposition.reverse()

    return decomposition


def linear_data_reconstruction(decomposition, G_operators):
    """Implementation of data reconstruction from the wavelet coeffs (phi_J, bar_phi_J, .., bar_phi_1).""" 
    # Reconstruct by a downward cascade starting with phi_J
    # following the red procedure in Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    phi = decomposition[0] 
    for i in range(len(G_operators)):
        bar_phi = decomposition[i+1]
        G_lo, G_hi = G_operators[i]
        phi = G_lo.conjugate(False).transpose(copy=False).dot(phi) + G_hi.conjugate(False).transpose(copy=False).dot(bar_phi)
        print(f'G_lo.H.shape {G_lo.conjugate(False).transpose(copy=False).shape}, bar_phi.shape {bar_phi.shape}, bar_phi count {bar_phi.shape[0]*bar_phi.shape[1]}, G_lo count_nonzero {np.count_nonzero(G_lo.toarray())}')
    # For clarity we highlight that final phi is on the lowest (finest) level
    # following Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    phi_0 = phi

    return phi_0
