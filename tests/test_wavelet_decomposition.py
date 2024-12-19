import numpy as np
import matplotlib.pyplot as plt
from hst.mrw1d import generate_wavelet, generate_G_operators, save_G_operators, data_decomposition, data_reconstruction

# Multi-coefficient super-Gaussian data generation
def superGaussian(x, cs): 
    data_length = len(x)
    Ndata = cs.shape[0]
    input_data = np.zeros((data_length, Ndata))
    for i in range(Ndata):
        input_data[:, i] = cs[i, 0]*np.exp(-((x-cs[i, 1])/cs[i, 2])**(2*cs[i, 3]))
    return input_data

# Save G_operators
save_Gs = False#True
# Print out essential properties of G_operators
verify_Gs = False#True

Nlevels = 7
Ndata = 15 # Number of randomized super-Gaussian (small number to plot)
data_length = 2**Nlevels * 100 # constructed as a proper  power of 2

#Nlevels = 3
#Ndata = 10
#data_length = 1024 # constructed as a proper  power of 2

# Generate data
domain_size = 1.0
x = np.linspace(-domain_size/2.0, domain_size/2.0, data_length)
# Random super-Gaussian c0*exp(-((x-c1)/c2)^(2*c3))
input_coeffs = np.random.rand(Ndata, 4)
# Rescale coeff values, so c0 in (0, 1), c1 in (-domain_size/2, domain_size/2), c2 in (0, domain_size/2), and c3 in {1, .., 10}
for i in range(Ndata):
    # c0 amplitude scaling is already correct
    input_coeffs[i, 0] = input_coeffs[i, 0]
    # c1 shift scaling
    input_coeffs[i, 1] = 0.8 * domain_size * (input_coeffs[i, 1] - 0.5)
    # c2 width scaling
    input_coeffs[i, 2] = 0.125 * domain_size * input_coeffs[i, 2]
    # c3 super-Gaussian power scaling
    input_coeffs[i, 3] = int(10 * input_coeffs[i, 3])
input_data = superGaussian(x, input_coeffs)
#input_data = np.random.rand(data_length, Ndata)

print(f'Nlevels {Nlevels}, data_length {data_length}, Ndata {Ndata}')

# Data decomposition into wavelet basis
wavelets = ['db1', 'db2', 'sdw2']
for wavelet in wavelets:
    print(f'TEST: wavelet {wavelet}')
    # Generate orthogonal G_lo (aka G) and G_hi (aka bar_G) operators
    G_operators = generate_G_operators(wavelet, Nlevels, data_length)
    print('G_operators generated.')
    if (save_Gs):
        # Save generated G_operators using dictionary
        file_name = f'G_operators_{wavelet}.npz'
        save_G_operators(G_operators, file_name)
        print(f'G_operators saved into {file_name}.')
    # Generate data decomposition into (phi_J, bar_phi_J, .., bar_phi_1)
    print('Compute data decomposition:')
    decomposition = data_decomposition(G_operators, input_data, verify_Gs)
    print('Decomposition done!')
    # Coarse phi_J interpretation
    cA = decomposition[0]
    #print(f'wavelet {wavelet}, cA')
    #print(f'{cA}')

    # Generate data reconstruction from (phi_J, bar_phi_J, .., bar_phi_1) 
    print('Compute data reconstruction:')
    reconstructed_data = data_reconstruction(decomposition, G_operators)
    print('Reconstruction done!')

    # Verification of direct and inverse multiresolution decomposition and reconstruction
    print(f'Verfication of the multiresolution decomposition and inverse reconstruction:')
    print(f'input_data.shape {input_data.shape}, reconstructed_data.shape {reconstructed_data.shape}')
    print('input_data[:, :] - reconstructed_data[:, :]')
    print(f'{input_data[:, :] - reconstructed_data[:, :]}')

    # Plot results
    for i in range(Ndata):
        plt.plot(x, input_data[:, i])
        plt.plot(x, reconstructed_data[:, i], '-.')
    plt.show()

# Discrete Wavelets analysis
wavelets = [] 
# DAUB4 or db2
dec_lo = 1.0 / 2.0**0.0 * np.array([-0.12940952255126037+0j, 0.2241438680420134+0j, 0.8365163037378079+0j, 0.48296291314453416+0j])
wavelets.append(['db2', dec_lo, [(4.0/3.0)**0.5]])

# Symmetric (a_k = a_{1-k}) SWD2 with indexes [-2, -1, 0, 1, 2, 3]
a1 = 0.662912 + 0.171163j
a2 = 0.110485 - 0.085581j
a3 = -0.066291 - 0.085581j

####################################################
# Orthogonal and invertibel scaling of BC for SDW2 #
####################################################
# Find the analytical solution for delta
from sympy import symbols, solve, lambdify, Matrix, conjugate

s_a1 = symbols('a_1')
s_a2 = symbols('a_2')
s_a3 = symbols('a_3')
s_alpha = symbols('alpha')
s_beta = symbols('beta')
s_gamma = symbols('gamma')
s_delta = symbols('delta')
G = Matrix([[s_alpha*s_a1, s_beta*s_a1, s_a2, s_a3, 0, 0, 0, 0], [s_a3, s_a2, s_a1, s_a1, s_a2, s_a3, 0, 0], [0, 0, s_a3, s_a2, s_a1, s_a1, s_a2, s_a3], [0, 0, 0, 0, s_a3, s_a2, s_beta*s_a1, s_alpha*s_a1]])
barG = Matrix([[-s_gamma*conjugate(s_a1), s_delta*conjugate(s_a1), -conjugate(s_a2), conjugate(s_a3), 0, 0, 0, 0], [-conjugate(s_a3), conjugate(s_a2), -conjugate(s_a1), conjugate(s_a1), -conjugate(s_a2), conjugate(s_a3), 0, 0], [0, 0, -conjugate(s_a3), conjugate(s_a2), -conjugate(s_a1), conjugate(s_a1), -conjugate(s_a2), conjugate(s_a3)], [0, 0, 0, 0, -conjugate(s_a3), conjugate(s_a2), -s_delta*conjugate(s_a1), s_gamma*conjugate(s_a1)]])

def print_corner(A, model):
    print(f'{model}')
    A00 = A[0, 0]
    A01 = A[0, 1]
    A10 = A[1, 0]
    A11 = A[1, 1]
    num_A00 = lambdify([s_a1, s_a2, s_a3], A00, 'numpy')
    num_A01 = lambdify([s_a1, s_a2, s_a3], A01, 'numpy')
    num_A10 = lambdify([s_a1, s_a2, s_a3], A10, 'numpy')
    num_A11 = lambdify([s_a1, s_a2, s_a3], A11, 'numpy')
    print(f'A[0, 0] {A00}')# = {num_A00(a1, a2, a3)}')
    print(f'A[0, 1] {A01}')# = {num_A01(a1, a2, a3)}')
    print(f'A[1, 0] {A10}')# = {num_A10(a1, a2, a3)}')
    print(f'A[1, 1] {A11}')# = {num_A11(a1, a2, a3)}')
    return A00, A01, A10, A11, num_A00, num_A01, num_A10, num_A11 

# Same s_alpha, s_beta multiplier of G and barG - DOES NOT WORK
#Og00, Og01, Og10, Og11, num_Og00, num_Og01, num_Og10, num_Og11 = print_corner(G * (barG.subs([(s_gamma, s_alpha), (s_delta, s_beta)])).H, 'G * barG.H')
#Iv00, Iv01, Iv10, Iv11, num_Iv00, num_Iv01, num_Iv10, num_Iv11 = print_corner(G.H * G + (barG.subs([(s_gamma, s_alpha), (s_delta, s_beta)])).H * barG.subs([(s_gamma, s_alpha), (s_delta, s_beta)]), 'G.H * G + barG.H * barG')

# s_alpha, s_beta multiplier of G and s_alpha^*, s_beta^* multiplier of barG - WORKS!!!
Og00, Og01, Og10, Og11, num_Og00, num_Og01, num_Og10, num_Og11 = print_corner(G * (barG.subs([(s_gamma, conjugate(s_alpha)), (s_delta, conjugate(s_beta))])).H, 'G * barG.H')
Iv00, Iv01, Iv10, Iv11, num_Iv00, num_Iv01, num_Iv10, num_Iv11 = print_corner(G.H * G + (barG.subs([(s_gamma, conjugate(s_alpha)), (s_delta, conjugate(s_beta))])).H * barG.subs([(s_gamma, conjugate(s_alpha)), (s_delta, conjugate(s_beta))]), 'G.H * G + barG.H * barG')

# Constrain is to satisfy orthogonality
equations = [Og00, Og01] # equal zero
solutions = solve(equations, s_alpha, s_beta, dict=True)
sol_number = 1 # positive solution #1 much better than negative solution #0
num_s_alpha = lambdify([s_a1, s_a2, s_a3], solutions[sol_number][s_alpha], 'numpy')
num_s_beta = lambdify([s_a1, s_a2, s_a3], solutions[sol_number][s_beta], 'numpy')
print(f's_alpha {solutions[sol_number][s_alpha]} = {num_s_alpha(a1, a2, a3)}')
print(f's_beta {solutions[sol_number][s_beta]} = {num_s_beta(a1, a2, a3)}')

alpha = num_s_alpha(a1, a2, a3)
beta = num_s_beta(a1, a2, a3)

print(f'alpha = {alpha}')
print(f'beta = {beta}')
####################################################

# Symmetry and unitary scaling 1/sqrt(2)
scale = 1.0 / 2.0**0.0
dec_lo = scale * np.array([a3, a2, a1, a1, a2, a3], dtype=complex)
wavelets.append(['SDW2', dec_lo, [alpha, beta]])

for wavelet, dec_lo, scaling in wavelets:
    dec_hi = generate_wavelet(dec_lo)
    
    print(f'{wavelet} filters')
    print(f'dec_lo {dec_lo}')
    print(f'dec_hi {dec_hi}')
    #print(f'Orthogonality')
    #print(f'dec_lo.dec_hi {dec_lo @ dec_hi}')
    #print('Invertibility')
    #print(f'dec_lo^*.dec_lo + dec_hi^*.dec_hi = {dec_lo.conjugate() @ dec_lo + dec_hi.conjugate() @ dec_hi}')

    print('HST') 
    Nrows = 4; Ncols = 2 * Nrows
    Npoints = len(dec_lo)
    mod = int(Npoints / 2)
    mat_G = np.zeros((Nrows, Ncols), dtype=complex)
    bar_mat_G = np.zeros((Nrows, Ncols), dtype=complex)
    for i in range(Nrows):
        for p in range(Npoints):
            #print(f'i {i}, 2*i+Npoints-mod-p {2*i+Npoints-mod-p}')
            # BC sets zero values of data outside the grid
            if (2*i+Npoints-mod-p >= 0 and 2*i+Npoints-mod-p <= Ncols-1):
                # Note that dec_lo and dec_hi are assigned in reverse order (compatible with Pywavelets)
                mat_G[i, 2*i+Npoints-mod-p] = dec_lo[p]
                bar_mat_G[i, 2*i+Npoints-mod-p] = dec_hi[p]

    # Apply the tailored scaling
    print(f'scaling of BC coefficients (first and last row): {scaling}')
    for j in range(len(scaling)):
        mat_G[0, j] = scaling[j] * mat_G[0, j]
        mat_G[Nrows-1][Ncols-1-j] = scaling[j] * mat_G[Nrows-1][Ncols-1-j] 
        bar_mat_G[0, j] = scaling[j].conjugate() * bar_mat_G[0, j]
        bar_mat_G[Nrows-1][Ncols-1-j] = scaling[j].conjugate() * bar_mat_G[Nrows-1][Ncols-1-j]

    mat_G = np.matrix(mat_G)
    bar_mat_G = np.matrix(bar_mat_G)
    print('Minimalistic mat_G with BCs')
    print(mat_G)
    print('Minimalistic bar_mat_G with BCs')
    print(bar_mat_G)
    print('Orthogonality mat_G.bar_mat_G^H')
    print(mat_G @ bar_mat_G.H)
    print('Invertibility mat_G^H.mat_G + bar_mat_G^H.bar_mat_G')
    print(mat_G.H @ mat_G + bar_mat_G.H @ bar_mat_G)
