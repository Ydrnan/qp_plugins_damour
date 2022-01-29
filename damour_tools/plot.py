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
    w_logx = False
    w_logy = False

    # states
    n_states = int(kwargs.get('n_states', 1))
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
        w_logx = True
       
    elif (kind == 'e=f(pt2)'):
        x = []
        y = []
        p_label = []
        for i in range(n_states):
            x.append(data[states[i]*5+2])
            y.append(data[states[i]*5+1])
            p_label.append(mo+"s st "+str(states[i]))
        x_label = 'PT2 energy (E$_h$)'
        y_label = 'Variational energy (E$_h$)'
        n_plots = n_states

    elif (kind == 'e=f(rpt2)'):
        x = []
        y = []
        p_label = []
        for i in range(n_states):
            x.append(data[states[i]*5+4])
            y.append(data[states[i]*5+1])
            p_label.append(mo+"s st "+str(states[i]))
        x_label = 'rPT2 energy (E$_h)$'
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
        w_logx = True
    else:
        print("Unknown kind of plot"+kind+" read the fuc*ing doc")
        exit()

    # Font
    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #rc('font',**{'family':'serif','serif':['Times']})
    #rc('text', usetex=True)
    for i in range(n_plots):
        if (kind == 'e=f(pt2)') or (kind == 'e=f(rpt2)'):
            ax.plot(x[i], y[i], label=p_label[i], marker='x', linewidth=2, markersize=6)
        else:
            ax.plot(x, y[i], label=p_label[i], marker='x', linewidth=2, markersize=6)

    ax.set_xlabel(x_label, fontsize=20)#, fontname='Times New Roman')
    #ax.set_xlim(left=0.1, right=2)
    ax.tick_params(axis='both',which='major',  labelsize=20, length=4, width=1)
    #ax.set_xticks(np.arange(0, 2.1, step=0.2))
    ax.set_ylabel(y_label, fontsize=20)#, fontname='Times New Roman')
    #ax.set_ylim(bottom=0.1, top=5.)
    if (w_logx): 
        ax.set_xscale('log')
    if (w_logy):
        ax.set_yscale('log')
    
    return ax

import sys

# test
import getopt

# init/ default
n_states = 1
basename = "plot"
title = None
grid = True
mticks = True

# read the command line arguments 
try:
    opts, args = getopt.getopt(sys.argv[1:],"hs:t:o:",["help","n_states=","title=","output=","grid=","mticks="])
except getopt.GetoptError:
    print("plot.py -h")
    exit()

for opt, arg in opts:
    if opt == "-h":
        print(" plot.py [options] <files>")
        print("         -s <Number of states>   , default=", n_states)
        print("         -t <Title>              , default=", title)
        print("         -o <Output name>        , default=", basename)
        print("         --grid=<Boolean>        , default=", grid)
        print("         --mticks=<Boolean>      , default=", mticks)
        sys.exit()
    #elif opt in ("-x", "--mol"):
    #    mol = arg
    #elif opt in ("-b", "--basis"):
    #    basis = arg
    elif opt in ("-s", "--n_states"):
        n_states = arg
        if isinstance(n_states,int):
            print("The number of states must be an integer")
            sys.exit()
            if n_states < 1:
                print("The number of states must be > 1")
                sys.exit()
    elif opt in ("-t", "--title"):
        title = arg
    elif opt in ("-o", "--output"):
        basename = arg
    elif opt in ("--grid"):
        grid = arg
        print(grid)
    elif opt in ("--mticks"):
        mticks = arg

# files
list_files = []

for arg in sys.argv:
    list_files.append(arg)
list_files.pop(0) # 1st time to remove the name of the prog

# remove options in list_files
i = 0
l = len(list_files)
while i < l:
    if list_files[i][:2] == '--':
        list_files.pop(i)
        l = l-1
    elif list_files[i][:1] == '-':
        list_files.pop(i)
        list_files.pop(i)
        l = l-2
    else:
        i = i+1
    print(i,l)

print(list_files)

list_methods = ['e=f(ndet)', 'e=f(pt2)', 'e=f(rpt2)', 'exc=f(ndet)']

for method in list_methods:
    figname = basename+'_'+method+'.png'
    print(figname)

    # Font
    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('font',**{'family':'serif','serif':['Times']})
    rc('text', usetex=True)

    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=False)

    for infile in list_files:
        print(infile)
        data = import_data(infile)
        print(n_states)
        print(grid)
        ax = my_plot(ax, data, kind=method, n_states=n_states, st_0=0, mo='OO', title=title)
 
    ax.legend(loc='best', fontsize=10)
    if (mticks):
        ax.minorticks_on()
    if (title != None):
        ax.set_title(title, fontsize=30)#, fontname='Times New Roman')
    if (grid):
        ax.grid(visible=True, which='major', axis='both', linewidth=0.5, linestyle='dashed')
    
    plt.subplots_adjust(left  = 0.16,  # The position of the left edge of the subplots, as a fraction of the figure width
        right = 0.98,   # The position of the right edge of the subplots, as a fraction of the figure width
        bottom = 0.125,  # The position of the bottom edge of the subplots, as a fraction of the figure height
        top = 0.9,      # The position of the top edge of the subplots, as a fraction of the figure height
        wspace = None,   # The width of the padding between subplots, as a fraction of the average Axes width
        hspace = None)   # The height of the padding between subplots, as a fraction of the average Axes width)
    plt.savefig(figname, dpi='figure', format='png', metadata=None,
            bbox_inches=None, pad_inches=0.1,
            backend=None)
    plt.show()


