import numpy as np
import matplotlib.pyplot as plt

# integral sqrt(a - b x^2)dx = 1/2 (x sqrt(a - b x^2) + (a b log(sqrt(a - b x^2) - sqrt(-b) x))/(-b)^(3/2)) + constant
#(assuming a complex-valued logarithm)

def W(q, E, m, omega):
    a = 2 * m * E
    b = m**2 * omega**2
    x = q
    return 1.0 / 2.0 * (x * (a - b * x**2)**0.5 + (a * b * np.log((a - b * x**2)**0.5 - (-b)**0.5 * x)) / (-b)**(3.0/2.0))

E = 1
m = 1
omega = 1

qmax = (2 * E / m)**0.5 / omega

t = np.linspace(0, 2.0 * np.pi / omega, 100)
q = qmax * np.sin(omega * t)

def S(t, q, E, m, omega):
    return W(q, E, m, omega) - E * t

S_evolution = S(t, q, E, m, omega)

plt.plot(t, S_evolution.real)
plt.plot(t, S_evolution.imag)
plt.show()
