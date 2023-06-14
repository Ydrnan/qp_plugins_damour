#!/bin/bash

source PATH_TO_QP2 #~/qp2/qp2_dev-stable/qp2/quantum_package.rc

ez=$1
frac=$2

if [ -z "$1" ]
then 
    echo "Missing 1st parameter: ezfio"
    exit
fi

if [ -z "$2" ]
then 
    echo "Missing 2nd parameter: fraction of correlation"
    exit
fi

re='^[0-9]+([.][0-9]+)?$'
if ! [[ $2 =~ $re ]] ; then
   echo "error: Not a number '$2'" >&2; exit 1
fi

if (( $(echo "$2 > 1.0" | bc -l) )) 
then
    echo "2nd parameter cannot be > 1"
    exit
fi

E_hf=$(cat ${ez}/hartree_fock/energy)
E_ccsd=$(cat ${ez}/ccsd/energy)
E_t=$(cat ${ez}/ccsd/energy_t)

percentage=$(perl -e "print $frac * 100")

if [ -z "$E_hf" ]
then
    echo "${ez}/hartree_fock/energy is empty"
    exit
fi

if [ -z "$E_ccsd" ]
then
    echo "${ez}/ccsd/energy is empty"
else
    target=$(perl -e "print $E_hf + $frac * ($E_ccsd - $E_hf)")
    echo ""
    echo "$percentage % of CCSD correlation: $target"
    echo "qp set percentage_exc target_energy_qmc $target"
fi

if [ -z "$E_t" ]
then
    echo "${ez}/ccsd/energy_t is empty"
else
    target=$(perl -e "print $E_hf + $frac * ($E_t - $E_hf)")
    echo ""
    echo "$percentage % of CCSD(T) correlation: $target"
    echo "qp set percentage_exc target_energy_qmc $target"
fi

