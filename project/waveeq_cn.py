import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as ct
import scipy.linalg as la
import time
import sim_params as par


def psi_init():
    """
    Metoda inicjująca funkcję falową
    :return: wartości początkowe funkcji falowej
    """
    psi = np.zeros(par.N + 1, dtype=complex)
    for ii in range(1, par.N):
        psi[ii] = np.exp(-(ii * par.dx - par.x0) ** 2 / (2 * par.sigma ** 2)) \
                  * np.exp(1.j * 2 * np.pi * (ii * par.dx - par.x0) / par.lam)
    norm_psi = psi / np.sqrt(norm(psi))
    return norm_psi


def v():
    """
    Metoda zwracająca wektor z wartościami potencjału
    :return: potencjał
    """
    _v = np.zeros(par.N + 1, dtype=complex)
    for ii in range(1, par.N):
        if 0.5 * par.N * par.dx < ii * par.dx < 0.5 * par.N * par.dx + par.l:
            _v[ii] = par.val*ct.e
    return _v


def norm(wave_fun):
    """
    Funkcja licząca normę funkcji falowej
    :param wave_fun: funkcja falowa
    :return: norma funkcji falowej
    """
    return np.trapz(np.abs(wave_fun)**2, dx=par.dx)


def kin_energy(wave_fun):
    """
    Funkcja licząca energie kinetyczną
    :param wave_fun: funkcja falowa
    :return: człon kinetyczny energii
    """
    tenergy = [0]
    for ii in range(1, len(wave_fun)-1):
        tenergy += [-ct.hbar**2/(2*ct.m_e)*(wave_fun[ii-1]-2*wave_fun[ii]+wave_fun[ii+1])/par.dx**2]
    tenergy += [0]
    return np.real(np.trapz(np.conj(wave_fun)*tenergy, dx=par.dx))/ct.e


def pot_energy(wave_fun, potential):
    """
    Funkcja licząca energię potencjalną
    :param wave_fun: funkcja falowa
    :param potential: potencjał
    :return: człon potencjalny energii
    """
    return np.real(np.trapz(np.abs(wave_fun) ** 2 * potential, dx=par.dx))/ct.e


def calculate_matrices(v):
    """
    Funkcja wyznaczająca macierze do użycia metody Crancka-Nicolsona
    :param v: lista wartości potencjału
    :return: a: macierz A z metody Crancka-Nicolsona
    :return: b: macierz B z metody Cranka-Nicolsona
    """

    alpha = -par.dt / (2 * par.dx * par.dx) * ct.hbar / (2 * ct.m_e * 1.j)
    a = np.zeros([par.N + 1, par.N + 1], dtype=complex)
    b = np.zeros([par.N + 1, par.N + 1], dtype=complex)

    a[0, 0] = a[par.N, par.N] = 1.
    b[0, 0] = b[par.N, par.N] = 1.
    for i in range(1, par.N):
        a[i, i] = 1 + 2 * alpha - 0.5* v[i]*par.dt/(1.j * ct.hbar)
        b[i, i] = 1 - 2 * alpha + 0.5* v[i]*par.dt/(1.j * ct.hbar)
        a[i, i + 1] = -alpha
        b[i, i + 1] = alpha
        a[i, i - 1] = -alpha
        b[i, i - 1] = alpha
    return a, b


def plot_func(xs, psi, tp, s_type, pot):
    """
    Funkcja do tworzenia wykresów
    :param xs: zdyskretyzowana oś OX
    :param psi: wartości funkcji w punktach z :param xs
    :param tp: wartości czasu
    :param s_type: typ symulacji (jedno-/wielowątkowa)
    :param pot: potencjał
    """
    plt.figure(s_type+': {0:.2f}'.format(tp*1e15)+" fs")
    plt.plot(xs, np.real(psi), xs, np.imag(psi), '-.')
    plt.legend(['$\Re{\psi}$', '$\Im{\psi}$'], loc='upper left')
    plt.xlabel("x")
    plt.ylabel("$\psi$")
    plt.axis([0, par.L, 1.5 * min(np.real(psi)), 1.5 * max(np.real(psi))])
    kine = kin_energy(psi)
    pote = pot_energy(psi, pot)
    plt.text(0, 0, "KE="+"{0:.4f}".format(kine))
    plt.text(35e-9, 0, "PE"+"{0:.4f}".format(pote))
    plt.text(20e-9, 0, "TE"+"{0:.4f}".format(pote+kine))
    plt.show()


def calc_psi_lu(wave_fun, lu_and_piv, b):
    """
    Funkcja do rozwiązywania równiania Schroedingera
    :param wave_fun: funkcja falowa
    :param lu_and_piv: macierze powstałe w wyniku rozkładu maczierzy A metodą LU
    :param b: macierz B z metody Crancka-Nicolsona
    :return:
    """
    temp_psi = la.lu_solve(lu_and_piv, np.dot(b, wave_fun))
    psi = temp_psi / np.sqrt(norm(temp_psi))
    return psi


def main_alg_lu(sim_type):
    """
    Funkcja główna z algorytmem przeprowadzającym symulację
    :param sim_type string z typem symulacji (wielo-/jednowątkowa)
    :return: całkowity czas rozwązywania równania metodą Cranka-Nicolsona
    """
    xs = np.arange(0, par.L + par.dx * 0.9, par.dx)
    t = 0.0
    psi = psi_init()
    tot_time = 0
    pot = v()
    a, b = calculate_matrices(pot)
    lu_decomp = la.lu_factor(a)
    while t < par.tmax:
        t1 = time.time()
        psi = calc_psi_lu(psi, lu_decomp, b)
        t2 = time.time()
        tot_time += t2 - t1
        for tp in par.tplot:
            if abs(t - tp) < 0.5 * par.dt:
                plot_func(xs, psi, tp, sim_type, pot)
        t += par.dt
    return tot_time
