import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from hst.wavelet_operators import generate_G_operators
from hst.linmultres import linear_data_decomposition, linear_data_reconstruction
from hst.nonlinmultres import nonlinear_data_decomposition, nonlinear_data_reconstruction
from hst.hst import hst_data_decomposition, hst_data_reconstruction
from hst.hst2 import hst2_data_decomposition, hst2_data_reconstruction
from hst.hst3 import hst3_data_decomposition, hst3_data_reconstruction
#from hst2 import hst2_data_decomposition, hst2_data_reconstruction

import matplotlib
def move_figure(f, x, y):
    """Move figure's upper left corner to pixel (x, y)"""
    backend = matplotlib.get_backend()
    if backend == 'tkagg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'wxagg':
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)
        
        
##############
# Input data #
##############
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

####################################
# General decomposition parameters #
####################################
# Symmetric Debauchy wavelet sdw2
wavelet = 'sdw2'

# Number of scatters
n_levels = 3 # = m + 1
print(f'Using wavelet {wavelet} on {n_levels} scattering levels within the Heisenberg scattering transform.')

# Generate orthogonal G_lo (aka G) and G_hi (aka bar_G) operators
G_operators = generate_G_operators(wavelet, n_levels, data_length)
print('G_operators generated.')

################################################################    
# linear Wavelet Transform decompostion and reconstruction #
################################################################
# Generate data decomposition into (S_J, bar_S_J, .., bar_S_1)
print('Compute data decomposition (S_J, bar_S_J, .., bar_S_1):')
decomposition_lwt = linear_data_decomposition(G_operators, input_data)
print('Decomposition done!')
# Coarse S_J interpretation
S_J_lwt = decomposition_lwt[0]
bar_S_J_lwt = decomposition_lwt[1]
bar_S_1_lwt = decomposition_lwt[n_levels]
print(f'S_J.shape {S_J_lwt.shape}, bar_S_J.shape {bar_S_J_lwt.shape}, bar_S_1.shape {bar_S_1_lwt.shape}')

# Generate data reconstruction from (S_J, bar_S_J, .., bar_S_1) 
print('Compute data reconstruction from (S_J, bar_S_J, .., bar_S_1):')
reconstructed_data_lwt = linear_data_reconstruction(decomposition_lwt, G_operators)
print('Reconstruction done!')

# Verification of direct and inverse multiresolution decomposition and reconstruction
print(f'Verfication of the nonlinear multiresolution decomposition and its inverse reconstruction:')
print(f'input_data.shape {input_data.shape}, reconstructed_data.shape {reconstructed_data_lwt.shape}')
print('input_data[:, :] - reconstructed_data[:, :]')
print(f'{input_data[:, :] - reconstructed_data_lwt[:, :]}')
print(f'norm: {np.linalg.norm(input_data[:, :] - reconstructed_data_lwt[:, :])}')


################################################################
# non-linear Wavelet Transform decompostion and reconstruction #
################################################################
# Unit operation on low-frequencies
def nonlinear_function(f):
    return f
def nonlinear_function_inverse(g):
    return g
# Logarithmic operation on high-frequencies
eps = 1e-10
c_nln = 1e-2
# Note that that the one shift is turned off
def R0(f):
    return f + c_nln #+ f / (abs(f) + eps)
def R0_inverse(g):
    return g - c_nln #- g / (abs(g) + eps)
def MG_bar_nonlinear_function(f):
    return 1j*np.log(R0(f))
def MG_bar_nonlinear_function_inverse(g): 
    return R0_inverse(np.exp(-1j*g))

# log scatter from S_0 to coarser structures
#def bar_nonlinear_function(f):
#    return MG_bar_nonlinear_function(f)
#def bar_nonlinear_function_inverse(g):
#    return MG_bar_nonlinear_function_inverse(g)

# log scatter from S_J to finer structures
def bar_nonlinear_function(f):
    return MG_bar_nonlinear_function(f) 
def bar_nonlinear_function_inverse(g):
    return MG_bar_nonlinear_function_inverse(g)

# Generate data decomposition into (S_J, bar_S_J, .., bar_S_1)
print('Compute data decomposition (S_J, bar_S_J, .., bar_S_1):')
decomposition_wt = nonlinear_data_decomposition(G_operators, input_data, nonlinear_function, bar_nonlinear_function)
print('Decomposition done!')
# Coarse S_J interpretation
S_J_wt = decomposition_wt[0]
bar_S_J_wt = decomposition_wt[1]
bar_S_1_wt = decomposition_wt[n_levels]
print(f'S_J.shape {S_J_wt.shape}, bar_S_J.shape {bar_S_J_wt.shape}, bar_S_1.shape {bar_S_1_wt.shape}')

# Generate data reconstruction from (S_J, bar_S_J, .., bar_S_1) 
print('Compute data reconstruction from (S_J, bar_S_J, .., bar_S_1):')
reconstructed_data_wt = nonlinear_data_reconstruction(decomposition_wt, G_operators, nonlinear_function_inverse, bar_nonlinear_function_inverse)
print('Reconstruction done!')

# Verification of direct and inverse multiresolution decomposition and reconstruction
print(f'Verfication of the nonlinear multiresolution decomposition and its inverse reconstruction:')
print(f'input_data.shape {input_data.shape}, reconstructed_data.shape {reconstructed_data_wt.shape}')
print('input_data[:, :] - reconstructed_data[:, :]')
print(f'{input_data[:, :] - reconstructed_data_wt[:, :]}')
print(f'norm: {np.linalg.norm(input_data[:, :] - reconstructed_data_wt[:, :])}')

#######################################
# HST decompostion and reconstruction #
#######################################
keep_bar_Sm = True

# Generate data decomposition into (S_J, bar_S_J, .., bar_S_1)
print('Compute data decomposition (S_J, bar_S_J, .., bar_S_1):')
decomposition = hst_data_decomposition(G_operators, input_data, keep_bar_Sm)
print('Decomposition done!')
# Coarse S_J interpretation
bar_S_m = decomposition[0]
S_m = decomposition[1]
S_0 = decomposition[n_levels - int(not keep_bar_Sm)]
print(f'S_J.shape {S_m.shape}, bar_S_J.shape {bar_S_m.shape}, S_0.shape {S_0.shape}')

# Generate data reconstruction from (S_J, bar_S_J, .., bar_S_1) 
print('Compute data reconstruction from (S_J, bar_S_J, .., bar_S_1):')
reconstructed_data = hst_data_reconstruction(decomposition, G_operators, keep_bar_Sm)
print('Reconstruction done!')

# Verification of direct and inverse multiresolution decomposition and reconstruction
print(f'Verfication of the nonlinear multiresolution decomposition and its inverse reconstruction:')
print(f'input_data.shape {input_data.shape}, reconstructed_data.shape {reconstructed_data.shape}')
print('input_data[:, :] - reconstructed_data[:, :]')
print(f'{input_data[:, :] - reconstructed_data[:, :]}')
print(f'norm: {np.linalg.norm(input_data[:, :] - reconstructed_data[:, :])}')

#################
# Visualization #
#################
def visualize(decomposition, reconstructed_data, fig_position, suptitle, Sj_text, S_J_text, file_name):
    # Prepare x vectors for each level by "bination"
    x_levels = []
    xlevel = x
    for i in range(n_levels):
        xlevel = xlevel[::2]
        print(f'xlevel.shape {xlevel.shape}')
        # Add x coordinates for bar_S_i
        x_levels.append(xlevel)
    # Reverse the order to make bar_S_J x coords first
    x_levels.reverse()

    # Plot results
    n_plots = 5
    dsnp = int(n_data / n_plots)
    print(f'dsnp {dsnp}')
    fig1, ax = plt.subplots()
    move_figure(fig1, fig_position[0], fig_position[1])
    fig2, axs = plt.subplots(n_plots, 3)
    move_figure(fig2, fig_position[0]+25, fig_position[1]+50)

    i_plot = 0
    for i in range(n_data):
        if i % dsnp == 0:
            print(f'Snapshot {i}')
            ax.plot(x, input_data[:, i].real)
            ax.plot(x, reconstructed_data[:, i].real, '-.')
            # S_J
            #ax_R = axs[i_plot, 0].twinx()
            
            if (keep_bar_Sm or Sj_text == "bar_Sj"):                
                axs[i_plot, 0].plot(x_levels[0], decomposition[0][:, i].real, 'x-', label='real')
                axs[i_plot, 0].plot(x_levels[0], decomposition[0][:, i].imag, 'o-', label='imag')
                
            # bar_S_j
            for j in range(n_levels):
                # Skip the S_J
                if (keep_bar_Sm or Sj_text == "bar_Sj"):
                    component = j + 1
                else:
                    component = j                    
                axs[i_plot, 1].plot(x_levels[j], decomposition[component][:, i].real, 'x-', label=f'j={n_levels-j}')
                axs[i_plot, 2].plot(x_levels[j], decomposition[component][:, i].imag, 'o-', label=f'j={n_levels-j}') 
            i_plot = i_plot + 1
    axs[0, 0].set_title(f'{S_J_text}, levels={n_levels}')
    axs[0, 0].legend()
    axs[0, 1].set_title(f'Re({Sj_text})')
    axs[0, 1].legend()
    axs[0, 2].set_title(f'Im({Sj_text})')
    axs[0, 2].legend()
    fig1.suptitle(suptitle)
    fig2.suptitle(suptitle)
    
#    plt.show()

os.makedirs('images', exist_ok=True)
visualize(decomposition_lwt, reconstructed_data_lwt, [100, 50], 'linear WT', "bar_Sj", "S_J", "images/lwt.png")
visualize(decomposition_wt, reconstructed_data_wt, [950, 50], 'non-linear WT', "bar_Sj", "S_J", "images/nwt.png")
visualize(decomposition, reconstructed_data, [1750, 50], 'HST', "Sj", "bar_S_J", "images/hst.png")


## HST2

decompositionHst2, decompositionHst2Full = hst2_data_decomposition(G_operators, input_data)
reconstructedHst2 = hst2_data_reconstruction(decompositionHst2, G_operators)

#################
# Visualization #
#################
full_tree_idx = 14

def visualize2(decomposition, decompositionFull, fig_position, reconstructed_data, suptitle, Sj_text, file_name):
    n_decomp = len(decomposition)
    
    # Prepare x vectors for each level by "bination"
    x_levels = []
    xlevel = x
    x_levels.append(xlevel)
    for i in range(n_levels):
        xlevel = xlevel[::2]
        #print(f'xlevel.shape {xlevel.shape}')
        # Add x coordinates for bar_S_i
        x_levels.append(xlevel)
    # Reverse the order to make bar_S_J x coords first
    #x_levels.reverse()

    # Plot results
    n_plots = 5
    dsnp = int(n_data / n_plots)
    print(f'dsnp {dsnp}')
    fig1, ax = plt.subplots() 
    move_figure(fig1, fig_position[0], fig_position[1])
    fig2, axs = plt.subplots(n_plots, 2)
    move_figure(fig2, fig_position[0]+25, fig_position[1]+50)

    i_plot = 0
    for i in range(n_data):
        if i % dsnp == 0 and i < input_data.shape[1]:
            print(f'Snapshot {i}')
            ax.plot(x, input_data[:, i].real)
            ax.plot(x, reconstructed_data[:, i].real, '-.')
            # S_J
            #ax_R = axs[i_plot, 0].twinx()
                            
            # bar_S_j
            for j in range(n_decomp):
                axs[i_plot, 0].plot(x_levels[-1], decomposition[j][:, i].real, 'x-', label=f'k={j}')
                axs[i_plot, 1].plot(x_levels[-1], decomposition[j][:, i].imag, 'o-', label=f'k={j}')
            i_plot = i_plot + 1
    axs[0, 0].set_title(f'Re({Sj_text})')
    axs[0, 0].legend()
    axs[0, 1].set_title(f'Im({Sj_text})')
    axs[0, 1].legend()
    fig1.suptitle(suptitle)
    fig1.tight_layout()
    fig2.suptitle(suptitle)
    fig2.tight_layout()
    
    fig3, axs2 = plt.subplots((n_levels+1), 2**n_levels)
    move_figure(fig3, fig_position[0]+50, fig_position[1]+100)
    for i in range(n_levels+1):
        for j in range(2**i):
            bar = '$\\bar{' if j % 2 == 1 else ''
            barend = '}$' if j % 2 == 1 else ''
            print(i, j, len(decompositionFull[i]))
            axs2[i, 2**(n_levels - i)*j].plot(x_levels[i], decompositionFull[i][j][:, full_tree_idx].real, 'x-', label=rf'Re({bar}S{barend}_{i}{j})')
            axs2[i, 2**(n_levels - i)*j].plot(x_levels[i], decompositionFull[i][j][:, full_tree_idx].imag, 'o-', label=rf'Im({bar}S{barend}_{i}{j})')
            axs2[i, 2**(n_levels - i)*j].legend()
    
    fig3.suptitle(f'{suptitle} full tree')
    fig3.tight_layout()
    
visualize2(decompositionHst2, decompositionHst2Full, [100, 700], reconstructedHst2, 'HST2', "Sj, bar_S_J", "images/hst2.png")


## HST3

decompositionHst3, decompositionHst3Full, _ = hst3_data_decomposition(G_operators, input_data)
reconstructedHst3 = hst3_data_reconstruction(decompositionHst3, G_operators)
    
visualize2(decompositionHst3, decompositionHst3Full, [950, 700], reconstructedHst3, 'HST3', "Sj, bar_S_J", "images/hst3.png")





def visualize3(decomposition, decompositionFull, decompositionFullWavelet, fig_position, reconstructed_data, suptitle, Sj_text, file_name):
    n_decomp = len(decomposition)
    
    # Prepare x vectors for each level by "bination"
    x_levels = []
    xlevel = x
    x_levels.append(xlevel)
    for i in range(n_levels):
        xlevel = xlevel[::2]
        #print(f'xlevel.shape {xlevel.shape}')
        # Add x coordinates for bar_S_i
        x_levels.append(xlevel)
    # Reverse the order to make bar_S_J x coords first
    #x_levels.reverse()

    # Plot results
    n_plots = 5
    dsnp = int(n_data / n_plots)
    print(f'dsnp {dsnp}')
    fig1, ax = plt.subplots() 
    move_figure(fig1, fig_position[0], fig_position[1])
    fig2, axs = plt.subplots(n_plots, 2)
    move_figure(fig2, fig_position[0]+25, fig_position[1]+50)

    i_plot = 0
    for i in range(n_data):
        if i % dsnp == 0 and i < input_data.shape[1]:
            print(f'Snapshot {i}')
            ax.plot(x, input_data[:, i].real)
            ax.plot(x, reconstructed_data[:, i].real, '-.')
            # S_J
            #ax_R = axs[i_plot, 0].twinx()
                            
            # bar_S_j
            for j in range(n_decomp):
                axs[i_plot, 0].plot(x_levels[-1], decomposition[j][:, i].real, 'x-', label=f'k={j}')
                axs[i_plot, 1].plot(x_levels[-1], decomposition[j][:, i].imag, 'o-', label=f'k={j}')
            i_plot = i_plot + 1
    axs[0, 0].set_title(f'Re({Sj_text})')
    axs[0, 0].legend()
    axs[0, 1].set_title(f'Im({Sj_text})')
    axs[0, 1].legend()
    fig1.suptitle(suptitle)
    fig1.tight_layout()
    fig2.suptitle(suptitle)
    fig2.tight_layout()
    
    fig3, axs2 = plt.subplots((n_levels+1), 2**n_levels)
    move_figure(fig3, fig_position[0]+50, fig_position[1]+100)
    for i in range(n_levels+1):
        for j in range(2**i):
            bar = '$\\bar{' if j % 2 == 1 else ''
            barend = '}$' if j % 2 == 1 else ''
            print(i, j, len(decompositionFull[i]))
            axs2[i, 2**(n_levels - i)*j].plot(x_levels[i], decompositionFull[i][j][:, full_tree_idx].real, 'x-', label=rf'Re({bar}S{barend}_{i}{j})')
            axs2[i, 2**(n_levels - i)*j].plot(x_levels[i], decompositionFull[i][j][:, full_tree_idx].imag, 'o-', label=rf'Im({bar}S{barend}_{i}{j})')
            axs2[i, 2**(n_levels - i)*j].plot(x_levels[i], decompositionFullWavelet[i][j][:, full_tree_idx].real, 'x--', label=rf'Re({bar}S{barend}_{i}{j})')
            axs2[i, 2**(n_levels - i)*j].plot(x_levels[i], decompositionFullWavelet[i][j][:, full_tree_idx].imag, 'o--', label=rf'Im({bar}S{barend}_{i}{j})')
            axs2[i, 2**(n_levels - i)*j].legend()
    
    fig3.suptitle(f'{suptitle} full tree')
    fig3.tight_layout()


## HST3 on other data

data_other = np.sin(np.linspace(-np.pi, np.pi, 1024))
input_data = [data_other, data_other]
input_data = np.array(input_data).T
full_tree_idx = 0

decompositionHst3, decompositionHst3Full, decompositionFullWavelet = hst3_data_decomposition(G_operators, input_data)
reconstructedHst3 = hst3_data_reconstruction(decompositionHst3, G_operators)
    
visualize3(decompositionHst3, decompositionHst3Full, decompositionFullWavelet, [1750, 700], reconstructedHst3, 'HST3 other', "Sj, bar_S_J", "images/hst3.png")

plt.show()