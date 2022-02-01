#!/bin/python3

import numpy as np

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
