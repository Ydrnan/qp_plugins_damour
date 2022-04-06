#!/bin/python3

import sys
from linear_regression import *
from convert_units import *
from import_data import *

def extrapolation_fci(filename, **kwargs):
    print("Filename:",filename)
    data = import_data(filename)
    
    n_states = np.size(data,0) // 5
    print("\nNumber of states:", n_states)

    n_data = np.size(data,1)

    if (n_data < 3):
        print("\nNeed at least 3 points, found "+str(n_data)+", exit")
        sys.exit()

    print("f(PT2) = a*PT2 + b\n")
    
    weight_type = kwargs.get('weight', 'None')
    print("Weight of the points for the linear regression:", weight_type)

    a = np.zeros((n_data, n_states))
    b = np.zeros((n_data, n_states))
    R2 = np.zeros((n_data, n_states))
    exc = np.zeros((n_data, n_states))

    # Extrapolated energies of the ground state for diff nb of points
    for state in range(0,n_states):
        E = data[1+5*state]
        PT2 = data[2+5*state]
        
        if weight_type == 'pt2': 
            weight = 1.0/PT2
        elif weight_type == 'pt2**2':
            weight = 1.0/PT2**2
        else:
            weight = np.full(len(E), 1.0)
    
        print("%s %2d" %("\nState:",state))
        if (state == 0):
            print("%3s %6s %13s %9s"%("n","a","b (Ha)","R^2"))
        else:
            print("%3s %6s %13s %12s %10s %7s"%("n","a","b (Ha)","ref (Ha)","Exc (eV)","R^2"))
    
        for nb_points in range(3,min(8,n_data)+1):
            i = nb_points-3
            # f(x) = a * x + b
            # with weights
            reg = lin_reg_v2(PT2,E,weight,nb_points)
            a[i][state] = reg[0] 
            b[i][state] = reg[1]
            R2[i][state] = reg[2]
    
            # without
            #reg = lin_reg(PT2,E,nb_points)
            #a=reg[0]
            #b=reg[1]
            #R2=reg[2]
    
            if (state == 0):
                print("| %1d | %6.4f | %8.4f | %8.6f |"%( nb_points, a[i][state], b[i][state], R2[i][state]))
            else:
                exc[i][state] = Ha_to_eV(b[i][state] - b[i][0])
                print("| %1d | %6.4f | %8.4f | %8.4f | %7.4f | %8.6f |"%( nb_points, a[i][state],
                      b[i][state], b[i][0], exc[i][state] ,R2[i][state]))

    # Energy with 4 points
    return n_data,a,b,R2,exc

if __name__ == "__main__":
    # for direct extrapolation from an input file
    import sys
    from linear_regression import *
    from convert_units import *
    from import_data import *
    
    # main for extrapolation
    import getopt
    
    # init/ default
    weight_type = 'none'
    
    # read the command line arguments 
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:w:",["help","file=","weight="])
    except getopt.GetoptError:
        print("python3 extrapolation_fci.py -h")
        sys.exit()
    
    for opt, arg in opts:
        if opt == "-h":
            print(" extrapolation_fci [options] <files>")
            print("         -f <File>          , default=", None)
            print("         --file=<File>      , default=", None)
            print("         -w <String>        , default=", weight_type)
            print("         --weight=<String>  , default=", weight_type)
            sys.exit()
    
        elif opt in ("-w", "--weight"):
            weight_type = arg
        elif opt in ("-f", "--file"):
            filename = arg
    
    print("Filename:",filename)
    
    res = extrapolation_fci(filename,weight=weight_type)
    print("\nFinished")

