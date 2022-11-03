#!/bin/python3

import os
import sys
import numpy as np
#from linear_regression import *

def lin_reg_v2(x,y,weight,nb_points):
    import numpy as np

    # last n points and reverse their order
    tmp_x = np.array(x[-nb_points:])
    tmp_x = tmp_x[::-1]

    tmp_y = np.array(y[-nb_points:])
    tmp_y = tmp_y[::-1]

    tmp_w = np.array(weight[-nb_points:])
    tmp_w = tmp_w[::-1]

    # linear regression
    fit = np.polynomial.polynomial.polyfit(tmp_x, tmp_y, deg=1, rcond=None, full=False, w=tmp_w)

    # f(x) = ax + b
    a = fit[1]
    b = fit[0]

    # bar_y = \sum_i^n tmp_y[i]/n
    bar_y = np.sum(tmp_y) / len(tmp_y)

    # hat_y[i] = a*tmp_x[i] + b 
    hat_y = a * tmp_x + b

    # SSR = \sum_i (tmp_y[i] - hat_y[i])**2
    SSR = np.sum((tmp_y-hat_y)**2)

    # SST = \sum_i (tmp_y[i] - bar_y)**2
    SST = np.sum((tmp_y - bar_y)**2)

    R2 = 1.0 - SSR/SST

    return(a,b,R2)

def poly_reg(x,y,weight,nb_points):
    import numpy as np

    # last n points and reverse their order
    tmp_x = np.array(x[-nb_points:])
    tmp_x = tmp_x[::-1]

    tmp_y = np.array(y[-nb_points:])
    tmp_y = tmp_y[::-1]

    tmp_w = np.array(weight[-nb_points:])
    tmp_w = tmp_w[::-1]

    # linear regression
    fit = np.polynomial.polynomial.polyfit(tmp_x, tmp_y, deg=2, rcond=None, full=False, w=tmp_w)

    # f(x) = ax^2 + bx + c
    a = fit[2]
    b = fit[1]
    c = fit[0]

    # bar_y = \sum_i^n tmp_y[i]/n
    bar_y = np.sum(tmp_y) / len(tmp_y)

    # hat_y[i] = a*tmp_x[i] + b 
    hat_y = a * tmp_x**2 + b * tmp_x + c

    # SSR = \sum_i (tmp_y[i] - hat_y[i])**2
    SSR = np.sum((tmp_y-hat_y)**2)

    # SST = \sum_i (tmp_y[i] - bar_y)**2
    SST = np.sum((tmp_y - bar_y)**2)

    R2 = 1.0 - SSR/SST

    return(a,b,c,R2)


def extract_dip(n_states,fname, **kwargs):

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
    remove_first = False
    if ndet[0] == 1 and n_states > 1:
        remove_first = True
        ndet = np.delete(ndet,0)

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
        if remove_first:
            for j in range(1,len(output)):
                dip[istate][j-1] = output[j].replace("\n","")
        else:
            for j in range(len(output)):
                dip[istate][j] = output[j].replace("\n","")

    # E 
    E = np.zeros((n_states,len(ndet)))
    for istate in range(0,n_states):
        stream = os.popen("grep '# E  ' "+fname+" | awk '{print $"+str(istate+3)+"}'")
        output = stream.readlines()
        if remove_first:
            for j in range(1,len(output)):
                E[istate][j-1] = output[j]
        else:
            for j in range(len(output)):
                E[istate][j] = output[j]

    # PT2
    pt2 = np.zeros((n_states,len(ndet)))
    for istate in range(0,n_states):
        stream = os.popen("grep '# PT2  ' "+fname+" | awk '{print $"+str(2*istate+3)+"}'")
        output = stream.readlines()
        if remove_first:
            for j in range(1,len(output)):
                pt2[istate][j-1] = output[j]
        else:
            for j in range(len(output)):
                pt2[istate][j] = output[j]

    # rPT2
    rpt2 = np.zeros((n_states,len(ndet)))
    for istate in range(0,n_states):
        stream = os.popen("grep '# rPT2  ' "+fname+" | awk '{print $"+str(2*istate+3)+"}'")
        output = stream.readlines()
        if remove_first:
            for j in range(1,len(output)):
                rpt2[istate][j-1] = output[j] 
        else:
            for j in range(len(output)):
                rpt2[istate][j] = output[j]
    
    # Out
    out_file = fname+"_w_dipoles.dat"
    f = open(out_file,"w")
    f.write('{:14}'.format("# Ndet"))
    f.write('{:16s}'.format("       E_0"))
    f.write('{:12s}'.format("  PT2_0"))
    f.write('{:12s}'.format("  rPT2_0"))
    f.write('{:12s}'.format("  ||µ_0||"))
    for istate in range(1,n_states):
        f.write('{:16s}'.format("       E_"+str(istate)))
        f.write('{:12s}'.format("  PT2_"+str(istate)))
        f.write('{:12s}'.format("  rPT2_"+str(istate)))
        f.write('{:12s}'.format("  Exc. "+str(istate)))
        f.write('{:12s}'.format("  f^l_"+str(istate)))
        f.write('{:12s}'.format("  f^v_"+str(istate)))
        f.write('{:12s}'.format("  f^m_"+str(istate)))
        f.write('{:12s}'.format(" ||µ_"+str(istate)+"||"))
    f.write("\n")

    for j in range(len(ndet)):
        f.write('{:12d}'.format(ndet[j]))
        f.write('{:16.8f}'.format(E[0][j]))
        f.write('{:12.8f}'.format(pt2[0][j]))
        f.write('{:12.8f}'.format(rpt2[0][j]))
        f.write('{:12.8f}'.format(dip[0][j]))
        for istate in range(1,n_states):
            f.write('{:16.8f}'.format(E[istate][j]))
            f.write('{:12.8f}'.format(pt2[istate][j]))
            f.write('{:12.8f}'.format(rpt2[istate][j]))
            f.write('{:12.8f}'.format(exc[istate-1][j]))
            f.write('{:12.8f}'.format(f_l[istate-1][j]))
            f.write('{:12.8f}'.format(f_v[istate-1][j]))
            f.write('{:12.8f}'.format(f_m[istate-1][j]))
            f.write('{:12.8f}'.format(dip[istate][j]))
        f.write("\n")
    f.close

    ############### Extrapolation ###################

    print("\nNumber of states:", n_states)

    n_data = len(ndet)

    if (n_data < 3):
        print("\nNeed at least 3 points, found "+str(n_data)+", exit")
        sys.exit()

    weight_type = kwargs.get('weight', 'None')
    print("Weight of the points for the linear/polynomial regression:", weight_type)

    ### Extrapolation µ = f(PT2) ###
    print("\n### Dipole moments ###")

    ##############################
    # Linear regression µ = f(PT2)
    print("\n# Linear regression, µ(PT2) = ax + b")
    a = np.zeros((n_data, n_states))
    b = np.zeros((n_data, n_states))
    R2 = np.zeros((n_data, n_states))
    exc = np.zeros((n_data, n_states))

     # Extrapolated dipole moments for diff nb of points
    for state in range(0,n_states):
        dip_i = dip[state][:] 
        PT2_i = pt2[state][:]

        if weight_type == 'pt2':
            weight = 1.0/PT2_i
        elif weight_type == 'pt2**2':
            weight = 1.0/PT2_i**2
        else:
            weight = np.full(len(ndet), 1.0)

        print("%s %2d %s" %("\nState:",state,", Dipole moment (Debye)"))
        print("%3s %10s %14s %9s"%("n","a","b (D)","R^2"))

        for nb_points in range(3,min(8,n_data)+1):
            i = nb_points-3
            # f(x) = a * x + b
            # with weights
            reg = lin_reg_v2(PT2_i,dip_i,weight,nb_points)
            a[i][state] = reg[0]
            b[i][state] = reg[1]
            R2[i][state] = reg[2]

            print("| %1d | %10.6f | %10.6f | %9.6f |"%( nb_points, a[i][state], b[i][state], R2[i][state]))
   
    ##################################
    # polynomial regression µ = f(PT2)
    print("\n# Polynomial regression, µ(PT2) = ax^2 + bx + c")
    a = np.zeros((n_data, n_states))
    b = np.zeros((n_data, n_states))
    c = np.zeros((n_data, n_states))
    R2 = np.zeros((n_data, n_states))
    exc = np.zeros((n_data, n_states))

     # Extrapolated dipole moments for diff nb of points
    for state in range(0,n_states):
        dip_i = dip[state][:]
        PT2_i = pt2[state][:]

        if weight_type == 'pt2':
            weight = 1.0/PT2_i
        elif weight_type == 'pt2**2':
            weight = 1.0/PT2_i**2
        else:
            weight = np.full(len(ndet), 1.0)

        print("%s %2d %s" %("\nState:",state,", Dipole moment (Debye)"))
        print("%3s %14s %14s %14s %9s"%("n","a","b","c (D)","R^2"))

        for nb_points in range(3,min(8,n_data)+1):
            i = nb_points-3
            # f(x) = a * x**2 + b * x + c
            # with weights
            reg = poly_reg(PT2_i,dip_i,weight,nb_points)
            a[i][state] = reg[0]
            b[i][state] = reg[1]
            c[i][state] = reg[2]
            R2[i][state] = reg[3]

            print("| %1d | %16.6f | %12.6f | %10.6f | %9.6f |"%( nb_points, a[i][state], b[i][state], c[i][state], R2[i][state]))

    ### Extrapolation oscillator strengths ###
    print("\n### Oscillator strengths ###")

    ##############################
    # Linear regression f = f(PT2)  
    print("\n# Linear regression, f(PT2) = ax + b")
    a = np.zeros((n_data, n_states))
    b = np.zeros((n_data, n_states))
    R2 = np.zeros((n_data, n_states))
    exc = np.zeros((n_data, n_states))

    # Extrapolated oscillator strengths for diff nb of points
    for gauge_type in range(0,3):
        for tr_state in range(0,n_states-1):

            PT2_i = (pt2[0][:] + pt2[tr_state+1][:])*0.5
            if gauge_type == 0:
                gauge = "length"
                f_i = f_l[tr_state][:]
            elif gauge_type == 1:
                gauge = "velocity"
                f_i = f_v[tr_state][:]
            else:
                gauge = "mixed"
                f_i = f_m[tr_state][:]
    
            if weight_type == 'pt2':
                weight = 1.0/PT2_i
            elif weight_type == 'pt2**2':
                weight = 1.0/PT2_i**2
            else:
                weight = np.full(len(ndet), 1.0)
   
            print("%s %2d %s %s %s" %("\nTransition n.",tr_state,", oscillator strength in",gauge,"gauge" ))
            print("%3s %10s %14s %9s"%("n","a","b  ","R^2"))

            for nb_points in range(3,min(8,n_data)+1):
                i = nb_points-3
                # f(x) = a * x + b
                # with weights
                reg = lin_reg_v2(PT2_i,f_i,weight,nb_points)
                a[i][tr_state] = reg[0]
                b[i][tr_state] = reg[1]
                R2[i][tr_state] = reg[2]
    
                print("| %1d | %10.6f | %10.6f | %9.6f |"%( nb_points, a[i][tr_state], b[i][tr_state], R2[i][tr_state]))

    ###################################
    # Polynomial regression, f = f(PT2)
    print("\n# Polynomial regression, f(PT2) = ax^2 + bx + c")
    a = np.zeros((n_data, n_states))
    b = np.zeros((n_data, n_states))
    c = np.zeros((n_data, n_states))
    R2 = np.zeros((n_data, n_states))
    exc = np.zeros((n_data, n_states))

    # Extrapolated oscillator strengths for diff nb of points
    for gauge_type in range(0,3):
        for tr_state in range(0,n_states-1):

            PT2_i = (pt2[0][:] + pt2[tr_state+1][:])*0.5
            if gauge_type == 0:
                gauge = "length"
                f_i = f_l[tr_state][:]
            elif gauge_type == 1:
                gauge = "velocity"
                f_i = f_v[tr_state][:]
            else:
                gauge = "mixed"
                f_i = f_m[tr_state][:]

            if weight_type == 'pt2':
                weight = 1.0/PT2_i
            elif weight_type == 'pt2**2':
                weight = 1.0/PT2_i**2
            else:
                weight = np.full(len(ndet), 1.0)

            print("%s %2d %s %s %s" %("\nTransition n.",tr_state,", oscillator strength in",gauge,"gauge" ))
            print("%3s %14s %14s %14s %9s"%("n","a","b","c  ","R^2"))

            for nb_points in range(3,min(8,n_data)+1):
                i = nb_points-3
                # f(x) = a * x**2 + b * x + c
                # with weights
                reg = poly_reg(PT2_i,f_i,weight,nb_points)
                a[i][tr_state] = reg[0]
                b[i][tr_state] = reg[1]
                c[i][tr_state] = reg[2]
                R2[i][tr_state] = reg[3]

                print("| %1d | %16.6f | %12.6f | %10.6f | %9.6f |"%( nb_points, a[i][tr_state], b[i][tr_state], c[i][tr_state], R2[i][tr_state]))


if __name__ == "__main__":

    fname = sys.argv[1]
    n_states = sys.argv[2]
    print(fname)
    print(n_states)
    res = extract_dip(n_states,fname,weight="none")
    print("\n=========================================================")
    print("#########################################################")
    print("=========================================================")
    res = extract_dip(n_states,fname,weight="pt2**2")
