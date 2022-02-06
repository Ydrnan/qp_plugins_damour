#!/bin/python3
from error_fci import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

# plot

# Font
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)

# init
fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=False)

epsilon = [1e-3,5e-4,1e-4,5e-5,1e-5]

for eps in epsilon:
    data = error_fci(eps)
    time = data[1]
    N_init = data[2]
    N = N_init
    time_T = time.T
    for j in range(0,3):
        N = N_init
        time = time_T[:][j]
        # remove the zero elements
        i = 0
        l = len(N)
        while i < l:
            if time[i] < 1e-12:
                time = np.delete(time,i)
                N = np.delete(N,i)
                l = l - 1
            else:
                i = i + 1
        print("new time", time)
        if (j == 0):
            ax.plot(N, time, label="$\epsilon^-=$"+str(eps),marker='x', linewidth=2, markersize=6)
        elif (j==2):
            ax.plot(N, time, label="$\epsilon^+=$"+str(eps),marker='x', linewidth=2, markersize=6)
        else:
            ax.plot(N, time, label="$\epsilon=$"+str(eps),marker='x', linewidth=2, markersize=6)

N = N_init[:7]
print(N)
ax.plot(N, np.exp(N), label="exp(N)",marker='x', linewidth=2, markersize=6)

ax.legend(loc='best', fontsize=10)
ax.set_yscale('log')
ax.set_xlabel("N(H)", fontsize=20)
ax.set_ylabel("Wall time (s)", fontsize=20)
print("\n epsilon", eps,"\nTime:\n",time)
plt.savefig("result_t_f_N_epsilon.png", dpi='figure', format='png', metadata=None,
                bbox_inches=None, pad_inches=0.1,
                backend=None)
plt.show()
