import project.waveeq_cn as we
import multiprocessing as mp
import matplotlib.pyplot as plt
import project.sim_params as par
import time


if __name__ == '__main__':
    par.menu()
    times1 = we.main_alg_lu("Jednowątkowa")
    p = mp.Pool(processes=mp.cpu_count())
    times2 = p.map(we.main_alg_lu, ["Wielowątkowa"])
    print(times2[0], times1)
    #plt.show()