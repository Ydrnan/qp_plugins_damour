#!/bin/python3

import os
import sys
import shutil

def diff(a,b):
    c = abs(a-b)
    return c

method_file = "debug_str.txt"
test_list = "debug_list.txt"

with open(method_file, 'r') as f:
    data = f.readlines()
    data.pop(0)
    data_split = []
    for line in data:
        data_split.append(line.split())
    method_name = []
    for i in range(len(data_split)):
        method_name.append(data_split[i][0])
    method_str = []
    for i in range(len(data_split)):
        method_str.append(data_split[i][1:])
    for i in range(len(method_str)):
        method_str[i] = " ".join(method_str[i][:])

#print(data_split)
#print(method_name)
#print(method_str)
#print("")

with open(test_list, 'r') as f:
    data = f.readlines()
    data.pop(0)
    data_split = []
    for line in data:
        data_split.append(line.split())
    list_ezfio = []
    for i in range(len(data_split)):
        list_ezfio.append(data_split[i][0])
    list_method = []
    for i in range(len(data_split)):
        list_method.append(data_split[i][1])
    list_val = []
    for i in range(len(data_split)):
        list_val.append(float(data_split[i][2]))
    list_threshold = []
    for i in range(len(data_split)):
        list_threshold.append(float(data_split[i][3]))


#print(list_ezfio)
#print(list_method)
#print(list_val)

error = []
nb_error = 0
for i in range(len(list_val)):

    # Copy 
    if os.path.isdir(list_ezfio[i]):
        shutil.rmtree(list_ezfio[i])
    
    ezfio = list_ezfio[i].replace(".","_save.")
    bash_cp = "cp -r "+ezfio+" "+list_ezfio[i]
    stream = os.popen(bash_cp)

    j = 0
    while list_method[i] != method_name[j]:
      j = j+1
    bash_str = method_str[j]
    bash_str = bash_str.replace("ezfio",list_ezfio[i])

    print("Test "+str(i)+": "+list_ezfio[i]+" "+list_method[i])

    stream = os.popen(bash_str)
    output = float(stream.readline())
     
    error.append(diff(output,list_val[i]))
    if error[i] > list_threshold[i]:
        print(" FAILED, error:",error[i])
        nb_error = nb_error + 1
    else:
        print(" OK    , error:",error[i])

    # rm
    shutil.rmtree(list_ezfio[i])

print("")
if nb_error > 0:
    print("Failed with:", nb_error, "errors")
else:
    print("Successful")
