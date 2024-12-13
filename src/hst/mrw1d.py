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
    # CSR construction data
    row = np.zeros(Nrows * Npoints - 2 * (mod-1)) 
    col = np.zeros(Nrows * Npoints - 2 * (mod-1))
    data_lo = np.zeros(Nrows * Npoints - 2 * (mod-1))
    data_hi = np.zeros(Nrows * Npoints - 2 * (mod-1))
    index = 0
    for i in range(Nrows):
        for p in range(Npoints):
            #print(f'i {i}, 2*i+Npoints-mod-p {2*i+Npoints-mod-p}')
            # BC sets zero values of data outside the grid
            if (2*i+Npoints-mod-p >= 0 and 2*i+Npoints-mod-p <= Ncols-1):
                #G_lo[i, 2*i+Npoints-mod-p] = dec_lo[p]
                #G_hi[i, 2*i+Npoints-mod-p] = dec_hi[p]
                row[index] = i
                col[index] = 2*i+Npoints-mod-p
                data_lo[index] = dec_lo[p]
                data_hi[index] = dec_hi[p]
                index = index + 1
    # Using scipy sparse matrix
    G_lo = csr_matrix((data_lo, (row, col)), shape = (Nrows, Ncols))
    G_hi = csr_matrix((data_hi, (row, col)), shape = (Nrows, Ncols))
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
    return G_lo, G_hi

def generate_G_operators(wavelet, Nlevels, data_length):
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

def verify_G_operators(G_operators):
    """Verify that the set of G_operators is orthogonal and invertible""" 
    for G_lo, G_hi in G_operators:
        print(f'Otrhogonality: G_lo.G_hi.T')
        print(f'{G_lo.dot(G_hi.T)}')
        print(f'Orthogonality: G_hi.G_lo.T')
        print(f'{G_hi.dot(G_lo.T)}')
        print(f'Invertibility: G_lo^T.G_lo + G_hi^T.G_hi = I')
        print(f'{G_lo.T.dot(G_lo) + G_hi.T.dot(G_hi)}')

def save_G_operators(G_operators, file_name):
    # Save generated G_operators using dictionary 
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

def data_decomposition(G_operators, data, verify_Gs=False):
    """Implementation of the wavelet decomposition (phi_J, bar_phi_J, .., bar_phi_1)"""
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

def data_reconstruction(decomposition, G_operators):
    """Implementation of data reconstruction from the wavelet coeffs (phi_J, bar_phi_J, .., bar_phi_1).""" 
    # Reconstruct by a downward cascade starting with phi_J
    # following the red procedure in Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    phi = decomposition[0]
    for i in range(len(G_operators)):
        bar_phi = decomposition[i+1]
        G_lo, G_hi = G_operators[i]
        phi = G_lo.T.dot(phi) + G_hi.T.dot(bar_phi)
        print(f'G_lo.T.shape {G_lo.T.shape}, bar_phi.shape {bar_phi.shape}, bar_phi count {bar_phi.shape[0]*bar_phi.shape[1]}, G_lo count_nonzero {np.count_nonzero(G_lo.toarray())}')
    # For clarity we highlight that final phi is on the lowest (finest) level
    # following Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    phi_0 = phi
    return phi_0
