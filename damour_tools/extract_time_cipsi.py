import os
import sys
from extract_E_cipsi import f_search_str

def extract_time_cipsi(fname):
    # out
    out_file = fname + ".tdat"
    
    # wall time after each cipsi iteration (davidson + pt2)
    stream = os.popen("grep -i 'Summary' -B 24 "+fname+" | grep 'WALL' | awk '{if(NR % 2)print $6}'")
    output = stream.readlines()
    
    time = []
    for elem in output:
        time.append(elem.replace("\n", ""))
    
    # call the function to search the string
    Ndet = f_search_str("Summary",fname)
    E    = f_search_str("# E ",fname)
    PT2  = f_search_str("# PT2",fname)
    rPT2 = f_search_str("# rPT2",fname)
    
    # write the output
    f = open(out_file,"w")
    # first line
    f.write('{:14s}'.format("    Ndet"))
    for i in range(0,1): #range(len(E[0])):
        f.write('{:3s}'.format(" E_"))
        f.write('{:12s}'.format(str(i)))
        f.write('{:12s}'.format('PT2'))
        f.write('{:12s}'.format('Error_PT2'))
        #f.write('{:12s}'.format('rPT2'))
        #f.write('{:12s}'.format('Error_rPT2'))
        f.write('{:12s}'.format('Time'))
    f.write("\n")
    # data
    for i in range(len(Ndet)):
        f.write('{:10d}'.format(int(Ndet[i][0])))
        for j in range(0,1): #len(E[0])):
            f.write('{:15f}'.format(float(E[i][j])))
            f.write('{:12f}'.format(float(PT2[i][2*j])))
            f.write('{:12f}'.format(float(PT2[i][2*j+1])))
            #f.write('{:12f}'.format(float(rPT2[i][2*j])))
            #f.write('{:12f}'.format(float(rPT2[i][2*j+1])))
            f.write('{:15f}'.format(float(time[i])))
        f.write("\n")
    f.close()

if __name__ == "__main__":

    fname = sys.argv[1]
    print(fname)
    res = extract_time_cipsi(fname)

