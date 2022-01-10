# for direct extrapolation from an input file
import sys
from linear_regression import *

filename = sys.argv[1]

print(filename)

# Ha to eV
def Ha_to_eV(x):
    return x*27.21138397128


f = open(filename,'r')
firstline = f.readline() # just some comments one the first line
lenght = len(firstline.split())
n_states = lenght//5
print("\n%17s %2d\n"%("Number of states:",n_states))
lines = f.readlines()
f.close

if (len(lines) < 3):
    print("Need at least 3 points, found "+str(len(lines))+", exit")
    exit()

print("f(PT2) = a*PT2 + b")
text = "n  a  b (E_h)   R^2"

ref = [] # Extrapolated energies of the ground state for diff nb of points
for state in range(0,n_states):
    E = []
    PT2 = []
    weight = []
    for line in lines:
        l = line.split()
        E.append(float(l[1+5*state]))
        PT2.append(float(l[2+5*state]))

    print("%s %2d" %("\nstate:",state))
    if (state == 0):
        print("%4s %10s %10s %14s"%("n","a","b (Ha)","R^2"))
    else:
        print("%4s %10s %10s %10s %10s %14s"%("n","a","b (Ha)","ref (Ha)","Exc (eV)","R^2"))

    i = 0
    for nb_points in range(3,min(8,len(lines))+1):
        for w in range(len(PT2)):
            weight.append(1/PT2[w])

        # f(x) = a * x + b
        # with weights
        #reg = lin_reg_v2(PT2,E,weight,nb_points)
        #a = reg[0] 
        #b = reg[1]
        #R2 = reg[2]

        # without
        reg = lin_reg(PT2,E,nb_points)
        a=reg[0]
        b=reg[1]
        R2=reg[2]

        if (state == 0):
            ref.append(b)
            print("| %2d | %8.4f | %8.4f | %12.6f |"%( nb_points, a, b, R2))
        else:
            exc = b - ref[i]
            exc = Ha_to_eV(exc)
            print("| %2d | %8.4f | %8.4f | %8.4f | %8.4f | %12.6f |"%( nb_points, a, b, ref[i], exc ,R2))
            i = i +1

