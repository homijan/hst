import pywt
from pywt import wavedec
import numpy as np
from scipy.sparse import csr_matrix
#from cupyx.scipy.sparse import csr_matrix

def one_level_G_operators(Nrows, dec_lo, dec_hi):
    """Generate orthogonal G operators with clamped BC given low and high filters"""
    Ncols = 2 * Nrows
    # Number of points in the filters
    Npoints = len(dec_lo)
    mod = int(Npoints / 2)
    G_lo = np.zeros((Nrows, Ncols))
    G_hi = np.zeros((Nrows, Ncols))
    for i in range(Nrows):
        for p in range(Npoints):
            #print(f'i {i}, 2*i+Npoints-mod-p {2*i+Npoints-mod-p}')
            # BC sets zero values of data outside the grid
            if (2*i+Npoints-mod-p >= 0 and 2*i+Npoints-mod-p <= Ncols-1):
                G_lo[i, 2*i+Npoints-mod-p] = dec_lo[p]
                G_hi[i, 2*i+Npoints-mod-p] = dec_hi[p]
    # Compensate for the clamped BC to ensure G*G^T + bar_G*bar_G^T = I
    if (mod == 2):
        # mod == 2 means a 4 point filter, where only three points are used 
        # because of the BC (one value was dropped), hence the scaling sqrt(4/3)
        scale = (4.0/3.0)**0.5
        i = 0; j = 0
        G_lo[i, j] = scale * G_lo[i, j]
        G_hi[i, j] = scale * G_hi[i, j]
        i = Nrows-1; j = Ncols-1
        G_lo[i, j] = scale * G_lo[i, j]
        G_hi[i, j] = scale * G_hi[i, j]
    # Using scipy sparse matrix
    G_lo = csr_matrix(G_lo)
    G_hi = csr_matrix(G_hi)
    return G_lo, G_hi

def G_operators(wavelet, Nlevels, data_length):
    """Generate full set of mutli-resolution G_operators (wavelet + binate)"""
    # Obtain low and high resolution wavelet filters
    dec_lo = pywt.Wavelet(wavelet).dec_lo
    dec_hi = pywt.Wavelet(wavelet).dec_hi

    G_operators = []
    Nrows = data_length
    for level in range(Nlevels):
        if (Nrows % 2 != 0):
            print(f'Cannot binate data of length {Nrows} at level {level}')
            quit()
        Nrows = int(Nrows / 2)
        G_lo, G_hi = one_level_G_operators(Nrows, dec_lo, dec_hi)
        G_operators.append([G_lo, G_hi]) 

    # Reverse the order of G_operators levels to go downward
    G_operators.reverse()

    return G_operators

def G_operators_verification(G_operators):
    """Verify that the set of G_operators is orthogonal and invertible""" 
    for G_lo, G_hi in G_operators:
        print(f'Otrhogonality: G_lo.G_hi.T')
        print(f'{G_lo.dot(G_hi.T)}')
        print(f'Orthogonality: G_hi.G_lo.T')
        print(f'{G_hi.dot(G_lo.T)}')
        print(f'Invertibility: G_lo^T.G_lo + G_hi^T.G_hi = I')
        print(f'{G_lo.T.dot(G_lo) + G_hi.T.dot(G_hi)}')

def data_decomposition(G_operators, data, verify_Gs=False):
    """Implementation of the wavelet decomposition (phi_J, bar_phi_J, .., bar_phi_1)"""
    if (verify_Gs):
        # Verify orthogonality and invertibility of G_operators at all levels
        G_operators_verification(G_operators)

    # Execute upward decompostion from fine to coarse
    upward_decomposition = []
    # The G_operators levels need to be reversed upward
    for G_lo, G_hi in reversed(G_operators):
        print(f'G_lo.shape {G_lo.shape}, data.shape {data.shape}')
        # Using scipy sparse matrix
        phi = G_lo.dot(data)
        bar_phi = G_hi.dot(data)
        #print(f'phi.shape {phi.shape}, bar_phi.shape {bar_phi.shape}')
        upward_decomposition.append([phi, bar_phi])
        print(f'data count {data.shape[0]*data.shape[1]}')
        print(f'G_lo count_nonzero {np.count_nonzero(G_lo.toarray())}')
        # current level data vector
        data = phi

    # Construct downward decomposition from coarse to fine
    # as vector (phi_J, bar_phi_J, bar_phi_J-1,.., bar_phi_1) as in
    decomposition = []
    phi_J = upward_decomposition[len(upward_decomposition)-1][0]
    decomposition.append(phi_J)
    for phi, bar_phi in reversed(upward_decomposition):
        #print(f'bar_phi.shape {bar_phi.shape}')
        decomposition.append(bar_phi)

    return decomposition

def data_reconstruction(decomposition, G_operators):
    """Implementation of data reconstruction from the wavelet coeffs (phi_J, bar_phi_J, .., bar_phi_1)."""
    # Downward cascade starting with data = phi_J
    phi_J = decomposition[0] 
    print(f'phi_J,shape {phi_J.shape}')
    # Reconstruct
    data = phi_J
    for i in range(len(G_operators)):
        bar_phi = decomposition[i+1]
        G_lo, G_hi = G_operators[i]
        print(f'data.shape {data.shape}, bar_phi.shape {bar_phi.shape}, G_lo.T.shape {G_lo.T.shape}, G_hi.T.shape {G_hi.T.shape}')
        data = G_lo.T.dot(data) + G_hi.T.dot(bar_phi)
    return data
