#!/bin/python3

import numpy as np
from import_data import *
from extrapolation_fci import *
import sys
from scipy import interpolate

def error_fci(epsilon_search):
    Ndet_search = np.zeros(((38-2)//2,3))
    time_search = np.zeros(((38-2)//2,3))
    N = np.zeros((38-2)//2)
    
    first = []
    last = []
    N_e_fci_array = []
    N_delta_e_fci = []

    for i in range(2,38,2):
        n = (i-2)//2
        N[n] = i
        filename = "H_"+str('{0:02d}'.format(i))+"_1.8.fci.out.tdat"
    
        data_extra = extrapolation_fci(filename)
        
        e_fci = data_extra[2] 
        n_e = data_extra[0] # size of e_fci
        # take at max the 4 first energies
        e_fci = e_fci[:min(len(e_fci),4)]
        if len(e_fci) >= 2:
            for e in e_fci:
                delta_e_fci = 2.0*max(abs(e-e_fci[1]))
        else:
                delta_e_fci = 0.0

        print("E_FCI:", e_fci, delta_e_fci)
        e_fci = e_fci[min(len(e_fci),2)-1] # 2nd element if len >= 2 
        e_fci_pm_delta = [e_fci - delta_e_fci, e_fci, e_fci + delta_e_fci]
        # array contening the different e_fci used
        N_e_fci_array.append(e_fci)
        N_delta_e_fci.append(delta_e_fci)
        print("e_fci",e_fci)

        data = import_data(filename)
        
        Ndet = data[0]
        E = data[1]
        PT2 = data[2]
        time = data[4]
    
        print("\nTime:",time)
        print("\nNdet:",Ndet)
        for j in range(0,3):
            e_fci = e_fci_pm_delta[j]

            # epsilon(Ndet,N) =(e_fci -  (E(Ndet,N) + PT2(Ndet,N))) / e_fci
            epsilon = (e_fci - (E + PT2)) / e_fci

            # debug
            #print("\nEpsilon:",epsilon)
            #last.append(epsilon[-1])
            #first.append(epsilon[0])
    
            # interpolation of Ndet(epsilon,N)
            f = interpolate.interp1d(epsilon, Ndet)
            print(Ndet_search)
            try:
                Ndet_search[n][j]= f(epsilon_search)
            except ValueError:
                Ndet_search[n][j] = 0.0
                print("Error in Ndet interpolation")
            print(epsilon_search, Ndet_search[n][j], i)
    
            # interpolation of time(Ndet,N) = time(epsilon,N)
            g = interpolate.interp1d(Ndet, time)
            try:
                time_search[n][j] = g(Ndet_search[n][j])
            except ValueError:
                time_search[n][j] = 0.0
                print("Error in time interpolation")
            print(Ndet_search[n][j], time_search[n][j], i)
    
    print("\n Results")
    print("%8s %14s %12s %14s %12s %14s %12s %7s"%("epsilon","Ndet","Time","Ndet","Time","Ndet","Time","N"))
    for i in range(2,38,2):
        n = (i-2)//2
        print("%8.4f %15.4f %12.4f %15.4f %12.4f %15.4f %12.4f %6d"%(epsilon_search,
            Ndet_search[n][0], time_search[n][0],
            Ndet_search[n][1], time_search[n][1],
            Ndet_search[n][2], time_search[n][2],i))

    print("\n E fci",N_e_fci_array)
    print("\n N   E_fci    delta_E_fci")
    for n in range(2,38,2):
        i = int((n-2)//2)
        print("%2d %10.6f %10.6f"%(n, N_e_fci_array[i], N_delta_e_fci[i]))

    return Ndet_search, time_search, N    
    #print(Ndet_search)
    #print(last)
    #print(first)
    #print(time_search)

# prog
if __name__ == "__main__":
    epsilon_search = float(sys.argv[1])

    res = error_fci(epsilon_search)
    print(res)


