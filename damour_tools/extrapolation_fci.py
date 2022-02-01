# for direct extrapolation from an input file
import sys
from linear_regression import *
from convert_units import *

def extrapolation_fci(filename, **kwargs):
    f = open(filename,'r')
    firstline = f.readline() # just some comments one the first line
    lenght = len(firstline.split())
    n_states = lenght//5
    print("\n%17s %2d\n"%("Number of states:",n_states))
    lines = f.readlines()
    f.close
    
    if (len(lines) < 3):
        print("Need at least 3 points, found "+str(len(lines))+", exit")
        sys.exit()
    
    print("f(PT2) = a*PT2 + b\n")
    
    weight_type = kwargs.get('weight', 'None')
    print("Weight of the points for the linear regression:", weight_type)

    ref = [] # Extrapolated energies of the ground state for diff nb of points
    for state in range(0,n_states):
        E = []
        PT2 = []
        weight = []
        for line in lines:
            l = line.split()
            E.append(float(l[1+5*state]))
            PT2.append(float(l[2+5*state]))

        print("%s %2d" %("\nState:",state))
        if (state == 0):
            print("%3s %6s %13s %9s"%("n","a","b (Ha)","R^2"))
        else:
            print("%3s %6s %13s %12s %10s %7s"%("n","a","b (Ha)","ref (Ha)","Exc (eV)","R^2"))
    
        i = 0
        for nb_points in range(3,min(8,len(lines))+1):
            for w in range(len(PT2)):
                if weight_type == 'pt2': 
                    weight.append(1.0/PT2[w])
                elif weight_type == 'pt2**2':
                    weight.append(1.0/PT2[w]**2)
                else:
                    weight.append(1.0)
    
            # f(x) = a * x + b
            # with weights
            reg = lin_reg_v2(PT2,E,weight,nb_points)
            a = reg[0] 
            b = reg[1]
            R2 = reg[2]
    
            # without
            #reg = lin_reg(PT2,E,nb_points)
            #a=reg[0]
            #b=reg[1]
            #R2=reg[2]
    
            if (state == 0):
                ref.append(b)
                print("| %1d | %6.4f | %8.4f | %8.6f |"%( nb_points, a, b, R2))
            else:
                exc = b - ref[i]
                exc = Ha_to_eV(exc)
                print("| %1d | %6.4f | %8.4f | %8.4f | %7.4f | %8.6f |"%( nb_points, a, b, ref[i], exc ,R2))
                i = i +1

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
