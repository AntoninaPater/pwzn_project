import matplotlib.pyplot as plt
import numpy as np
from math import pi
import multiprocessing as mp
import time

m = 9.109e-31
hbar = 1.0545718e-34
L = 40e-9
dx = 5e-9
N = int(L / dx)
tmax = 100e-15
dt = 0.02e-15
sigma = 25 * dx
ev = 1.6e-19

def psi_init(n, x0=10e-9, lam=5e-9):
    """
    :param lam: długość fali -> p = hbar * 2 * PI / lambda
    :param n: liczba węzłów
    :param x0: początkowa wwspółrzędna położenia
    :return: wartości początkowe funkcji falowej
    """
    psi = np.zeros(n + 1, dtype=complex)
    print("energia", (hbar*2*np.pi/lam)**2/2/m/ev)
    for ii in range(1, n):
        psi[ii] = np.exp(-(ii * dx - x0) ** 2 / (2 * sigma ** 2)) * np.exp(1.j * 2 * np.pi * (ii * dx - x0) / lam)
    norm_psi = psi / np.sqrt(norm(psi))
    print("kin", kin_energy(norm_psi))
    return norm_psi


def v(n, l, val):
    """
    :param val: wartość potencjału w eV
    :param l: długość progu
    :param n: liczba węzłów

    :return: potencjał
    """
    _v = np.zeros(n + 1, dtype=complex)
    for ii in range(1, n):
        if 0.5 * n * dx < ii * dx < 0.5 * n * dx + l:
            _v[ii] = val*ev
    return _v


def norm(wave_fun):
    return np.trapz(np.abs(wave_fun)**2, dx=dx)


def kin_energy(wave_fun):
    tenergy = [0]
    for ii in range(1, len(wave_fun)-1):
        tenergy += [-hbar**2/(2*m)*(wave_fun[ii-1]-2*wave_fun[ii]+wave_fun[ii+1])/dx**2]
    tenergy += [0]
    return np.real(np.trapz(np.conj(wave_fun)*tenergy, dx=dx))/ev


def pot_energy(wave_fun, potential):
    return np.real(np.trapz(np.abs(wave_fun) ** 2 * potential, dx=dx))/ev


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
        b[i, i] = 1 - 2 * alpha + v[i]*dt/(1.j * hbar)
        a[i, i + 1] = -alpha
        b[i, i + 1] = alpha
        a[i, i - 1] = -alpha
        b[i, i - 1] = alpha
    return a, b


def plot_func(xs, psi, tp):
    """

    :param xs: węzły
    :param psi: wartości funkcji w więzłach
    :param tp: wartości czasu
    """
    plt.figure()
    plt.plot(xs, np.real(psi), xs, np.imag(psi), '-.', label="t = " + str(tp) + " s")
    plt.xlabel("x")
    plt.ylabel("$\psi$")
    plt.axis([0, 40e-9, 1.5 * min(np.real(psi)), 1.5 * max(np.real(psi))])
    plt.legend(loc="upper left")


def calc_psi(wave_fun, a, b):
    temp_psi = np.linalg.solve(a, np.dot(b, wave_fun))
    psi = temp_psi / np.sqrt(norm(temp_psi))
    return psi


def main_fun():
    tplot = [0, 90e-15]
    xs = np.arange(0, L + dx * 0.9, dx)
    t = 0.0
    psi = psi_init(N, lam=3e-9)
    tot_time = 0
    pot = v(N, L / 2, 0.1)
    a, b = calculate_matrices(pot)
    while t < tmax:
        t1 = time.time()
        psi = calc_psi(psi, a, b)
        t2 = time.time()
        tot_time += t2 - t1
        for tp in tplot:
            if abs(t - tp) < 0.5 * dt:
                plot_func(xs, psi, tp)
                kine = kin_energy(psi)
                pote = pot_energy(psi, pot)
                print("KE: ", kine, '\n', "PE: ", pote, '\n', "TE: ", kine + pote, sep='')
        t += dt
    print(tot_time)
    plt.show()


if __name__ == '__main__':
    p = mp.Pool(mp.cpu_count())
    p.map(main_fun, [])
    print("OK")
    main_fun()
