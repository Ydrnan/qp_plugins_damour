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
    f_l=np.zeros((n_states-1,len(ndet)))
    f_v=np.zeros((n_states-1,len(ndet)))
    f_m=np.zeros((n_states-1,len(ndet)))
    exc=np.zeros((n_states-1,len(ndet)))

    for istate in range(1,n_states+1):
        j = 0
        for line in read_file:
            if line.startswith("   #  Transition n.{:3d}:".format(istate)):
                line = line.split()
                f_l[istate-1][j] = float(line[11].replace(",",""))
                f_v[istate-1][j] = float(line[13].replace(",",""))
                f_m[istate-1][j] = float(line[15].replace(",",""))
                exc[istate-1][j] = float(line[5])
                j = j + 1
        load_file.close()

    # Dipole moment
    dip=np.zeros((n_states,len(ndet)))
    for istate in range(0,n_states):
        stream = os.popen("grep -A "+str(n_states+1)+" 'Dipole moments (D)' "+fname+" | grep ' "+str(istate)+" ' | awk '{print $5}'")
        output = stream.readlines()
        for j in range(len(output)):
            dip[j][istate] = output[j].replace("\n","")
    
    # Out
    out_file = fname+"_w_dipoles.dat"
    f = open(out_file,"w")
    f.write('{:12s}'.format("# Ndet"))
    f.write('{:10s}'.format("||µ_0||"))
    for istate in range(1,n_states):
        f.write('{:10s}'.format("Exc. "+str(istate)))
        f.write('{:10s}'.format("f^l_"+str(istate)))
        f.write('{:10s}'.format("f^v_"+str(istate)))
        f.write('{:10s}'.format("f^m_"+str(istate)))
        f.write('{:10s}'.format("||µ_"+str(istate)+"||"))
    f.write("\n")

    for j in range(len(ndet)):
        f.write('{:10d}'.format(ndet[j]))
        f.write('{:10.6f}'.format(dip[0][j]))
        for istate in range(1,n_states):
            f.write('{:10.6f}'.format(exc[istate-1][j]))
            f.write('{:10.6f}'.format(f_l[istate-1][j]))
            f.write('{:10.6f}'.format(f_v[istate-1][j]))
            f.write('{:10.6f}'.format(f_m[istate-1][j]))
            f.write('{:10.6f}'.format(dip[istate][j]))
        f.write("\n")
    f.close
    


if __name__ == "__main__":

    fname = sys.argv[1]
    n_states = sys.argv[2]
    print(fname)
    print(n_states)
    res = extract_dip(n_states,fname)
