import pywt
from pywt import wavedec
import numpy as np
from scipy.sparse import csr_matrix
#from cupyx.scipy.sparse import csr_matrix

def generate_wavelet(dec_lo):
    """Create a consistent high frequency decomposition wavelet b_k = (-1)^k a^*_{1-k}"""
    N = len(dec_lo)
    dec_hi = np.zeros(N, dtype=complex)
    for index in range(N):
        # offeset of the local k from python i index as k_loc = i + offset
        offset = int(1 - N / 2)
        k_loc = index + offset
        #print(f'array index {index}, k_loc {k_loc}, 1 - k_loc {1 - k_loc}')
        # b_k = (-1)^k a^*_{1-k}
        dec_hi[index] = (-1)**k_loc * dec_lo[1 - k_loc - offset].conjugate()
    return dec_hi

def get_wavelet(wavelet_type):
    if wavelet_type=='db1':
        dec_lo = pywt.Wavelet('db1').dec_lo
        scaling = []
    elif wavelet_type=='db2':
        dec_lo = pywt.Wavelet('db2').dec_lo
        scaling = [(4.0/3.0)**0.5]
    elif wavelet_type=='sdw2':
        # Symmetric (a_k = a_{1-k}) SWD2 with indexes [-2, -1, 0, 1, 2, 3]
        a1 = 0.662912 + 0.171163j
        a2 = 0.110485 - 0.085581j
        a3 = -0.066291 - 0.085581j
        dec_lo = np.array([a3, a2, a1, a1, a2, a3], dtype=complex)
        scaling = [1.0204969605748353+0.015876976569495486j, 1.012298185699968-0.015876906424713344j]
    dec_hi = generate_wavelet(dec_lo)
    # Assign dict
    wavelet = {}
    wavelet['mother'] = dec_lo
    wavelet['daughter'] = dec_hi
    wavelet['scaleBC'] = scaling
    return wavelet

def one_level_G_operators(Nrows, dec_lo, dec_hi, scaling):
    """Generate orthogonal G operators with clamped BC given low and high filters"""

    Ncols = 2 * Nrows
    # Number of points in the filters
    Npoints = len(dec_lo)
    mod = int(Npoints / 2)
    # CSR construction data
    row = np.zeros(Nrows * Npoints - 2 * (mod-1), dtype=complex) 
    col = np.zeros(Nrows * Npoints - 2 * (mod-1), dtype=complex)
    data_lo = np.zeros(Nrows * Npoints - 2 * (mod-1), dtype=complex)
    data_hi = np.zeros(Nrows * Npoints - 2 * (mod-1), dtype=complex)
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
    #print(f'scaling of BC coefficients (first and last row): {scaling}')
    for j in range(len(scaling)):
        G_lo[0, j] = scaling[j] * G_lo[0, j]
        G_lo[Nrows-1, Ncols-1-j] = scaling[j] * G_lo[Nrows-1, Ncols-1-j] 
        G_hi[0, j] = scaling[j].conjugate() * G_hi[0, j]
        G_hi[Nrows-1, Ncols-1-j] = scaling[j].conjugate() * G_hi[Nrows-1, Ncols-1-j]
    return G_lo, G_hi


def generate_G_operators(wavelet_type, Nlevels, data_length):
    """Generate full set of mutli-resolution G_operators (wavelet + binate)"""
    # Obtain low and high resolution wavelet filters and boundary scaling
    wavelet = get_wavelet(wavelet_type)
    dec_lo = wavelet['mother']
    dec_hi = wavelet['daughter']
    scaling = wavelet['scaleBC']
  
    G_operators = []
    Nrows = data_length
    for level in range(Nlevels):
        if (Nrows % 2 != 0):
            print(f'Cannot binate data of length {Nrows} at level {level}')
            quit()
        Nrows = int(Nrows / 2)
        G_lo, G_hi = one_level_G_operators(Nrows, dec_lo, dec_hi, scaling)
        G_operators.append([G_lo, G_hi]) 

    # Reverse the order of G_operators levels to go downward
    G_operators.reverse()

    return G_operators


def verify_G_operators(G_operators):
    """Verify that the set of G_operators is orthogonal and invertible""" 
    for G_lo, G_hi in G_operators:
        print(f'Otrhogonality: G_lo.G_hi.H')
        print(f'{G_lo.dot(G_hi.H)}')
        print(f'Orthogonality: G_hi.G_lo.H')
        print(f'{G_hi.dot(G_lo.H)}')
        print(f'Invertibility: G_lo^H.G_lo + G_hi^H.G_hi = I')
        print(f'{G_lo.H.dot(G_lo) + G_hi.H.dot(G_hi)}')


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
        phi = G_lo.H.dot(phi) + G_hi.H.dot(bar_phi)
        print(f'G_lo.H.shape {G_lo.H.shape}, bar_phi.shape {bar_phi.shape}, bar_phi count {bar_phi.shape[0]*bar_phi.shape[1]}, G_lo count_nonzero {np.count_nonzero(G_lo.toarray())}')
    # For clarity we highlight that final phi is on the lowest (finest) level
    # following Fig. 2 in Marchand et al, Wavelet Conditional Renormalization Group (2022)
    phi_0 = phi
    return phi_0
