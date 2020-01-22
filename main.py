import waveeq_cn
import multiprocessing as mp
import matplotlib.pyplot as plt
import sim_params as par
import time


if __name__ == '__main__':
    par.menu()
    times1 = waveeq_cn.main_alg_lu("Jednowątkowa")
    p = mp.Pool(processes=mp.cpu_count())
    times2 = p.map(waveeq_cn.main_alg_lu, ["Wielowątkowa"])
    print(times2[0], times1)
    #plt.show()