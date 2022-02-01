#!/bin/python3

# Prog to generate plots:

# files=$(ls -R */*fci.out)
# python extract_E_cispi $files
# files=$(ls -R */*.dat)
# python plot.py -s <Number of states> $files

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import matplotlib.font_manager

def parse_boolean(s):
    return s == 'True'

def Ha_to_eV(e):
    return e*27.211382543519

def import_data(filename, **kwargs):
    last = int(kwargs.get('last', 0))
    with open(filename,'r') as f:
        trash = f.readline() # comments line
        lines = f.readlines()
    
    data = []
    for line in lines:
        data.append(line.split())
    data = np.array(data)
    if (last != 0):
        n=len(data)
        data=data[n-last:][:]
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
    print("States: ",states)

    # checks the states ordering
    for i in range(1,n_states):
        if (states[i-1] >= states[i]):
            print("st_"+str(i-1)+" must be smaller than st_"+str(i)+": "+str(states[i-1])+">="+str(states[i])+", exit")
            exit()

    # kind
    if (kind == 'e=f(ndet)'):
        x = []
        y = []
        p_label = []
        for i in range(n_states):
            x.append(data[0])
            y.append(data[states[i]*5+1])
            p_label.append(mo+" st "+str(states[i]))
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
            p_label.append(mo+" st "+str(states[i]))
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
            p_label.append(mo+" st "+str(states[i]))
        x_label = 'rPT2 energy (E$_h)$'
        y_label = 'Variational energy (E$_h$)'
        n_plots = n_states

    elif (kind == 'exc=f(ndet)'):
        x = []
        y = []
        p_label = []
        if (n_states >= 2):
            for i in range(n_states-1):
                x.append(data[0])
                y.append(Ha_to_eV(data[states[i+1]*5+1]+data[states[i+1]*5+2]-data[states[i]*5+1]-data[states[i]*5+2]))
                p_label.append(mo+"s st "+str(states[i])+" -$>$ "+str(states[i+1]))
        else:
            print("n_states="+str(n_states)+", no excitation energy, exit")
            exit()
        x_label = 'Number of determinants'
        y_label = 'Excitation energy $\Delta_{E+PT2}$ (eV)'
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
        ax.plot(x[i], y[i], label=p_label[i], marker='x', linewidth=2, markersize=6)

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
title = ""
grid = True
mticks = True
show=False

# read the command line arguments 
try:
    opts, args = getopt.getopt(sys.argv[1:],"hs:t:o:",["help","n_states=","title=","output=","grid=","mticks=","show="])
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
        print("         --show=<Boolean>        , default=", show)
        sys.exit()
    #elif opt in ("-x", "--mol"):
    #    mol = arg
    #elif opt in ("-b", "--basis"):
    #    basis = arg
    elif opt in ("-s", "--n_states"):
        try:
            n_states = int(arg)
        except ValueError: 
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
        try:
            grid = parse_boolean(arg)            
        except ValueError:
            print("The argument for the grib must be a boolean")
            sys.exit()
    elif opt in ("--mticks"):
        try:
            mticks = parse_boolean(arg)
        except ValueError:
            print("The argument for the mticks must be a boolean")
            sys.exit()
    elif opt in ("--show"):
        try:
            show = parse_boolean(arg)
        except ValueError:
            print("The argument for show must be a boolean")
            sys.exit()

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

print("Files:", list_files)
print("N states:", n_states)

if n_states >= 2:
    list_methods = ['e=f(ndet)', 'e=f(pt2)', 'e=f(rpt2)', 'exc=f(ndet)']
else:
    list_methods = ['e=f(ndet)', 'e=f(pt2)', 'e=f(rpt2)']

# 0 -> takes all points, n /= 0 -> takes the n last points
zoom = [0,7]

# loop over the methods/kinds
for method in list_methods:

    # loop to zoom
    for last in zoom:
        if (last == 0):
            figname = basename+'_'+method+'.png'
        else:
            figname = basename+'_'+method+'_'+'zoom'+str(last)+'.png'

        # Font
        #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        rc('font',**{'family':'serif','serif':['Times']})
        rc('text', usetex=True)

        # init
        fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=False)

        # loop over the files to add the plots
        for infile in list_files:
            mo = infile.replace(".fci.out.dat","")
            mo = mo.replace(".opt_fci.out.dat","")
            mo = mo.replace("_"," ")
            title = title.replace("_"," ")
            data = import_data(infile, last=last)
            ax = my_plot(ax, data, kind=method, n_states=n_states, st_0=0, st_1=1, mo=mo, title=title)
 
        ax.legend(loc='best', fontsize=10)
        if (mticks):
            ax.minorticks_on()
        if (title != None):
            ax.set_title(title, fontsize=30)#, fontname='Times New Roman')
        ax.grid(visible=grid, which='major', axis='both', linewidth=0.5, linestyle='dashed')
        # margins
        plt.subplots_adjust(left  = 0.17,  # The position of the left edge of the subplots, as a fraction of the figure width
            right = 0.98,   # The position of the right edge of the subplots, as a fraction of the figure width
            bottom = 0.125,  # The position of the bottom edge of the subplots, as a fraction of the figure height
            top = 0.9,      # The position of the top edge of the subplots, as a fraction of the figure height
            wspace = None,   # The width of the padding between subplots, as a fraction of the average Axes width
            hspace = None)   # The height of the padding between subplots, as a fraction of the average Axes width)
        # Save
        plt.savefig(figname, dpi='figure', format='png', metadata=None,
                bbox_inches=None, pad_inches=0.1,
                backend=None)
        if (show):
            plt.show()


