#!/bin/python3

import os
import sys
import numpy as np

def extract_dip(n_states,fname):

    n_states = int(n_states)
    ndet = []
    load_file = open(fname,"r")
    read_file = load_file.readlines()

    # Det
    for line in read_file:
        if line.startswith("Summary"):
            line = line.split()
            
            ndet.append(int(line[4]))
    load_file.close()

    ndet = np.array(ndet)

    # Oscillator strength and excitation energy
    f_l=np.zeros((len(ndet),n_states-1))
    f_v=np.zeros((len(ndet),n_states-1))
    f_m=np.zeros((len(ndet),n_states-1))
    exc=np.zeros((len(ndet),n_states-1))
    for i in range(1,n_states+1):
        j = 0
        for line in read_file:
            if line.startswith("   #  Transition n.{:3d}:".format(i)):
                line = line.split()
                f_l[j][i-1] = float(line[11].replace(",",""))
                f_v[j][i-1] = float(line[13].replace(",",""))
                f_m[j][i-1] = float(line[15].replace(",",""))
                exc[j][i-1] = float(line[5])
                j = j + 1
        load_file.close()

    # Dipole moment
    dip=np.zeros((len(ndet),n_states))
    for istate in range(0,n_states):
        stream = os.popen("grep -A "+str(n_states+1)+" 'Dipole moments (D)' "+fname+" | grep ' "+str(istate)+" ' | awk '{print $5}'")
        output = stream.readlines()
        for i in range(len(output)):
            dip[i][istate] = output[i].replace("\n","")
    
    # Out
    out_file = fname+"_w_dipoles.dat"
    f = open(out_file,"w")
    f.write('{:12s}'.format("# Ndet"))
    f.write('{:10s}'.format("||µ_0||"))
    for istate in range(0,n_states-1):
        f.write('{:10s}'.format("Exc. "+str(istate)))
        f.write('{:10s}'.format("f^l_"+str(istate)))
        f.write('{:10s}'.format("f^v_"+str(istate)))
        f.write('{:10s}'.format("f^m_"+str(istate)))
        f.write('{:10s}'.format("||µ_"+str(istate)+"||"))
    f.write("\n")

    for i in range(len(ndet)):
        f.write('{:10d}'.format(ndet[i]))
        f.write('{:10.6f}'.format(dip[i][0]))
        for istate in range(0,n_states-1):
            f.write('{:10.6f}'.format(exc[i][istate]))
            f.write('{:10.6f}'.format(f_l[i][istate]))
            f.write('{:10.6f}'.format(f_v[i][istate]))
            f.write('{:10.6f}'.format(f_m[i][istate]))
            f.write('{:10.6f}'.format(dip[i][istate]))
        f.write("\n")
    f.close
    


if __name__ == "__main__":

    fname = sys.argv[1]
    n_states = sys.argv[2]
    print(fname)
    print(n_states)
    res = extract_dip(n_states,fname)
