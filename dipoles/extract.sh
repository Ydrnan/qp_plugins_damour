#!/bin/bash

#source ~/App/qp2/quantum_package.rc 

# ./extract.sh file n_states

file=$1
n_states=$2

# Ndet
echo "# Ndet      " > tmp_ndet.txt
grep "Summary at N_det" $file | awk '{printf "# %-10s\n", $5}' >> tmp_ndet.txt

# Extract dipole moment
state=0
while [ $state -lt $n_states ]
do
  echo "Dip. st. $state" > tmp_dip_${state}.txt
  grep -A $(($n_states+1)) "Dipole moments (D)" $file | grep " $state " | awk '{printf "%-.6f\n", $5}' >> tmp_dip_${state}.txt

  state=$(($state+1))
done

# Extract Exc.
tr=1
while [ $tr -lt $n_states ]
do
  echo "Exc. $tr(eV)" > tmp_exc_${tr}.txt
  grep "#  Transition n.  $tr" $file | awk '{printf "%-.6f\n", $6}' >>  tmp_exc_${tr}.txt

  tr=$(($tr+1))
done

# Extract oscillator strength 
tr=1
while [ $tr -lt $n_states ]
do
  echo "Osc. str. $tr" > tmp_o_${tr}.txt
  grep "#  Transition n.  $tr" $file | awk '{printf "%-.6f\n", $12}' >> tmp_o_${tr}.txt

  tr=$(($tr+1))
done

# Paste everything
paste tmp_ndet.txt tmp_dip_0.txt > tmp.txt
rm tmp_ndet.txt tmp_dip_0.txt

tr=1
while [ $tr -lt $n_states ]
do
  paste tmp.txt tmp_exc_${tr}.txt > tmp1_${tr}.txt
  rm tmp.txt tmp_exc_${tr}.txt
  paste tmp1_${tr}.txt tmp_o_${tr}.txt > tmp2_${tr}.txt
  rm tmp1_${tr}.txt tmp_o_${tr}.txt
  paste tmp2_${tr}.txt tmp_dip_${tr}.txt > tmp.txt
  rm tmp2_${tr}.txt tmp_dip_${tr}.txt
  tr=$(($tr+1))
done

# Extract Ndet, E, PT2
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/extract_E_cipsi.py -f $file

# Add data in comments
echo "" >> ${file}.dat
cat tmp.txt >> ${file}.dat

# Plot file
sed -i "s/#//g" tmp.txt
sed -i "2d" tmp.txt
mv tmp.txt ${file}_w_dipoles.dat

