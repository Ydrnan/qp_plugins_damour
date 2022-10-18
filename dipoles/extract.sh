#!/bin/bash

#source ~/App/qp2/quantum_package.rc 

# ./extract.sh file n_states

file=$1
n_states=$2

# Extract Ndet, E, PT2
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/extract_E_cipsi.py -f $file
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/extrapolation_fci.py -f $file.dat > $file.extrapolation_fci.dat
python3 $QP_ROOT/plugins/qp_plugins_damour/dipoles/extract_dip.py $file $n_states > $file.extrapolation_dipole.dat

cp ${file}_w_dipoles.dat tmp.txt
sed -i 's/^/#/' tmp.txt

# Add data in comments
echo "" >> ${file}.dat
cat tmp.txt >> ${file}.dat

rm tmp.txt
