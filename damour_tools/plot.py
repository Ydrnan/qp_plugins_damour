#!/bin/python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import matplotlib.font_manager

infile = 'cn3_loc_fb_opt.opt_fci.out.dat'


def import_data(filename):
    with open(filename,'r') as f:
        trash = f.readline() # comments line
        lines = f.readlines()
    
    data = []
    for line in lines:
        data.append(line.split())
    data = np.array(data)
    data = data.T
    data = np.asarray(data, dtype=float, order='F') # str to float
    
    return data

def my_plot(ax, data, kind, **kwargs):
    # states
    n_states = kwargs.get('n_states', 1)
    mo = kwargs.get('mo', '')
    states=[]
    for i in range(n_states):
        str_st = 'st_'+str(i)
        states.append(kwargs.get(str_st, i))
        print(str_st+": "+str(states[i]))

    # checks the states ordering
    for i in range(1,n_states):
        if (states[i-1] >= states[i]):
            print("st_"+str(i-1)+" must be smaller than st_"+str(i)+": "+str(states[i-1])+">="+str(states[i])+", exit")
            exit()

    print(states)

    # kind
    if (kind == 'e=f(ndet)'):
        x = data[0]
        y = []
        p_label = []
        for i in range(n_states):
            y.append(data[states[i]*5+1])
            p_label.append(mo+"s st "+str(states[i]))
        x_label = 'Number of determinants'
        y_label = 'Variational energy (E$_h$)'
        n_plots = n_states
       
    elif (kind == 'e=f(pt2)'):
        x = data[2]
        y = []
        p_label = []
        for i in range(n_states):
            y.append(data[states[i]*5+2])
            p_label.append(mo+"s st "+str(states[i]))
        x_label = 'PT2 energy (E$_h$)'
        y_label = 'Variational energy (E$_h$)'
        n_plots = n_states
    elif (kind == 'e=f(rpt2)'):
        x = data[2]
        y = []
        p_label = []
        for i in range(n_states):
            y.append(data[states[i]*5+4])
            p_label.append(mo+"s st "+str(states[i]))
        x_label = 'rPT2 energy E$_h$'
        y_label = 'Variational energy (E$_h$)'
        n_plots = n_states
    elif (kind == 'exc=f(ndet)'):
        x = data[0]
        y = []
        p_label = []
        if (n_states >= 2):
            for i in range(n_states-1):
                y.append(data[states[i+1]*5+1]+data[states[i+1]*5+2]-data[states[i]*5+1]-data[states[i]*5+2])
                p_label.append(mo+"s st "+str(states[i])+" -$>$ "+str(states[i+1]))
        else:
            print("n_states="+str(n_states)+", no excitation energy, exit")
            exit()
        x_label = 'Number of determinants'
        y_label = 'Excitation energy $\Delta_{E+PT2}$ (E$_h)$'
        n_plots = n_states-1
    else:
        print("Unknown kind of plot"+kind+" read the fuc*ing doc")
        exit()

    # Font-
    #matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('font',**{'family':'serif','serif':['Times']})
    rc('text', usetex=True)
    
    # output
    fname = 'toto.png'
    
    # E = f(Ndet)
    #fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=False)
    for i in range(n_plots):
        ax.plot(x, y[i], label=p_label[i], marker='x', linewidth=2, markersize=6)

    ax.legend(loc='best', fontsize=10)
    ax.set_xlabel(x_label, fontsize=20)#, fontname='Times New Roman')
    #ax.set_xlim(left=0.1, right=2)
    ax.tick_params(axis='both',which='major',  labelsize=20, length=4, width=1)
    #ax.set_xticks(np.arange(0, 2.1, step=0.2))
    ax.set_ylabel(y_label, fontsize=20)#, fontname='Times New Roman')
    #ax.set_ylim(bottom=0.1, top=5.)
    #ax.set_yscale('log')
    ax.set_xscale('log')
    ax.minorticks_on()
    ax.set_title("A beautiful title", fontsize=30)#, fontname='Times New Roman')
    ax.grid(visible=True, which='major', axis='both', linewidth=0.5, linestyle='dashed')
    #    plt.subplots_adjust(left  = 0.16,  # The position of the left edge of the subplots, as a fraction of the figure width
    #        right = 0.98,   # The position of the right edge of the subplots, as a fraction of the figure width
    #        bottom = 0.125,  # The position of the bottom edge of the subplots, as a fraction of the figure height
    #        top = 0.9,      # The position of the top edge of the subplots, as a fraction of the figure height
    #        wspace = None,   # The width of the padding between subplots, as a fraction of the average Axes width
    #        hspace = None)   # The height of the padding between subplots, as a fraction of the average Axes width)
    #    plt.savefig(fname, dpi='figure', format='png', metadata=None,
    #            bbox_inches=None, pad_inches=0.1,
    #            backend=None)
    #    plt.show()
    #
    #data = import_data(infile)
    #
    #out = my_plot(data,kind='exc=f(ndet)', n_states=2, mo='OO')
    return ax

import sys

# files
list_files = []

for arg in sys.argv:
    list_files.append(arg)
list_files.pop(0)

print(list_files)

fname = 'test.png'
fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=False)

for infile in list_files:
    print(infile)
    data = import_data(infile)
    ax = my_plot(ax, data, kind='e=f(ndet)', n_states=1, st_0=0, mo='OO')
    #ax = my_plot(ax, data, kind='e=f(ndet)', n_states=1, st_0=1, mo='NO')

plt.subplots_adjust(left  = 0.16,  # The position of the left edge of the subplots, as a fraction of the figure width
    right = 0.98,   # The position of the right edge of the subplots, as a fraction of the figure width
    bottom = 0.125,  # The position of the bottom edge of the subplots, as a fraction of the figure height
    top = 0.9,      # The position of the top edge of the subplots, as a fraction of the figure height
    wspace = None,   # The width of the padding between subplots, as a fraction of the average Axes width
    hspace = None)   # The height of the padding between subplots, as a fraction of the average Axes width)
plt.savefig(fname, dpi='figure', format='png', metadata=None,
        bbox_inches=None, pad_inches=0.1,
        backend=None)
plt.show()


