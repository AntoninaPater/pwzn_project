import scipy.constants as ct
L = 40e-9
dx = 1e-10
N = int(L / dx)
tmax = 150e-15
dt = 0.02e-15
sigma = 2.5e-9
tplot = [0, 90e-15]
x0 = 10e-9
lam = 5e-9
val = 0.15
l = 20e-9


def check(strin, temp):
    """
    Funkcja sprawdzająca poprawność danych
    :param strin: string z nazwą zmiennej
    :param temp: sprawdzana wartość
    :return: True, jeśli poprawna, False wpp.
    """
    if strin == 'x0' and 0 <= temp <= 20:
        return True
    elif strin == 'lam' and temp >= 0:
        return True
    elif strin == 'sigma' and temp > 0:
        return True
    elif strin == 'val' and abs(val) < 5:
        return True
    elif strin == 'l' and 0 < temp < 20:
        return True
    else:
        return False


def selector(string, cond, coeff):
    """
    Funkcja do zapisywania parametrów symulacji
    :param string: komunikat dla użytkownika
    :param cond: zmienna, którą chcemy wpisać
    :param coeff: współczynnik zamiany jednostek
    :return:
    """
    while True:
        temp = input(string)
        if temp == '':
            print("Wartość domyślna")
            return eval(cond)
            break
        elif check(cond, float(temp)):
            var = float(temp) * coeff
            print("Wartość poprawna")
            return var
            break
        else:
            print("Wartość nieprawidłowa")


def time_selector():
    """
    Procedura służąca wyborowi wartości czasów, dla których rysowane są wykresy
    """
    print("Domyślnie funkcja falowa jest wykreślania dla czasu początkowego i czasu 90 fs.")
    global tplot
    while True:
        add = input("Ile dodatkowych czasów chcesz wprowadzić [0]?\n")
        try:
            if add == '':
                print("Start symulacji")
            else:
                ii = 0
                while ii < int(add):
                    try:
                        addt = float(input("Podaj czas w fs (od 0 do 100) \n"))
                        if 0 < addt < 100:
                            tplot += [addt * 1e-15]
                            ii += 1
                        else:
                            print("Nieprawidłowa wartość")
                    except ValueError:
                        print("Nieprawidłowa wartość")
                print("Start symulacji")
                tplot.sort()
            break
        except ValueError:
            print("Nieprawidłowa wartość")


def menu():
    """
    Procedura realizująca menu dla użytkownika
    """
    instr = input("Wykonać symulację z parametrami domyślnymi (tak/nie)? [nie]\n")
    if not instr.lower() == 'tak':
        global x0, sigma, lam, val, l
        x0 = selector("Podaj początkowe położenie cząstki w nm (od 0 do 20) [10]\n", 'x0', 1e-9)
        lam = selector("Podaj początkową długość fali związanej z cząstką w nm (większa od 0) [5]\n", 'lam', 1e-9)
        print("Przybliżona energia  kinetyczna cząstki to:", "{0:.4f}".format(ct.hbar**2*(2*ct.pi/lam)**2/(2*ct.m_e)/ct.e), "eV")
        sigma = selector("Podaj początkową szerokość paczki związanej z cząstką w nm (większa od 0) [2.5]\n",
                         'sigma', 1e-9)
        val = selector("Podaj wartość potencjału w eV (moduł mniejszy od 1) [0.15]\n", 'val', 1)
        l = selector("Podaj długość progu w nm (od 0 do 20) [20]", 'l', 1e-9)
        time_selector()




