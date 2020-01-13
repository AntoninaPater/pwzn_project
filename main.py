import matplotlib.pyplot as plt
import numpy as np
from math import pi

m = 9.109e-31
hbar = 1.0545718e-34
L = 40e-9
dx = 0.1e-9
N = int(L / dx)
print(N)
tmax = 70e-15
dt = 0.02e-15
sigma = 5e-9
k = 2e10


def psi_init(n, x0=15e-9, lam=5e-9):
    """
    :param lam: ?
    :param n: liczba węzłów
    :param x0: początkowa wwspółrzędna położenia
    :return: wartości początkowe funkcji falowej
    """
    psi = np.zeros(n + 1, dtype=complex)
    for ii in range(1, n):
        psi[ii] = np.exp(-(ii * dx - x0) ** 2 / (2 * sigma ** 2)) * np.exp(1.j * 2 * np.pi * (ii * dx - x0) / lam)
    return psi


def v(n, l, val):
    """
    :param val: wartość potencjału w eV
    :param l: długość progu
    :param n: liczba węzłów
    :param x: współrzędna x miejsca, w którym wyliczamy potencjał

    :return: potencjał
    """
    _v = np.zeros(n + 1, dtype=complex)
    for ii in range(1, n):
        if 0.5 * l < ii * dx < l:
            _v[ii] = val*1.6*10**-19
    return _v


def calculate_matrices(v):
    """
    :param v: lista wartości potencjału
    :return: A: macierz
    :return: B: macierz
    """
    alpha = -dt / (2 * dx * dx) * hbar / (2 * m * 1.j)
    a = np.zeros([N + 1, N + 1], dtype=complex)
    b = np.zeros([N + 1, N + 1], dtype=complex)

    a[0, 0] = a[N, N] = 1.
    b[0, 0] = b[N, N] = 1.
    for i in range(1, N):
        a[i, i] = 1 + 2 * alpha
        b[i, i] = 1 - 2 * alpha +(v[i]*dt)/(1.j*pi)
        a[i, i + 1] = -alpha
        b[i, i + 1] = alpha
        a[i, i - 1] = -alpha
        b[i, i - 1] = alpha
    return a, b


def plot_func(xs, psi, tp):
    """

    :param xs: więzły
    :param psi: wartości funkcji w więzłach
    :param tp: wartości czasu
    """
    plt.figure()
    plt.plot(xs, np.real(psi), xs, np.imag(psi), '--', label="t = " + str(tp) + " s")
    plt.xlabel("x")
    plt.ylabel("$\psi$")
    plt.axis([0, 40e-9, 1.5 * min(np.real(psi)), 1.5 * max(np.real(psi))])
    plt.legend(loc="upper left")
    plt.show()

    # def main_fun(): # to na dole do środka i ustawić parametry te z góry
if __name__ == '__main__':

    tplot = [0, 34e-15, 68e-15]
    xs = np.arange(0, L + dx * 0.9, dx)
    t = 0.0
    psi = psi_init(N)

    a,b = calculate_matrices(v(N, L/2, 1))
    while t < tmax:
        psi = np.linalg.solve(a, np.dot(b, psi))
        for tp in tplot:
            if abs(t - tp) < 0.5 * dt:
                plot_func(xs, psi, tp)
        t += dt

