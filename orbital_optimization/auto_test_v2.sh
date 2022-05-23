#!/bin/bash

source $QP_ROOT/quantum_package.rc

# script to check the results and compare them with the results of reference

if [ -d "reference_results" ] 
then
    rm -r reference_results
fi

tar zxf reference_results.tar.gz
cd reference_results
cp ../debug_test.py .
cp ../debug_str.txt .
cp ../debug_list.txt .

python3 debug_test.py

cd ..
rm -r reference_results
