#!/bin/bash
#SBATCH -p xeonv6 -N 1 -n 1 -c 32 --exclusive

set -e

export OMP_PROC_BIND=false
# Put the right path to source qp2
#source $QP_ROOT/quantum_package.rc
source ~/App/qp2/quantum_package.rc

#module list
echo "Hostename" $HOSTNAME
echo "OMP_PROC_BIND" $OMP_PROC_BIND
#module load intel/2019.3
#module load gcc/9.2.0
#module load python/3.7.6-gcc-9.2.0

# All the echo are just here to find easily where are the bug if
# the script doesn't work

XYZ=h2
MOL=${XYZ}_double
BASIS=cc-pvqz  
CHARGE=0
MULTIPLICITY=1

# Studied states
STATE_0=1
STATE_1=2

# To select the states if they are close 
SELECT_STATES=true       # if the states are close
N_MAX_STATES=4           # In this case, how many states ?
N_DET_MAX_SELECT=1e3     # And with how many determinants ?
N_DET_MAX_HF=1e3
N_DET_MAX_NO=1e3
N_DET_MAX_OO=1e3

# Method to break the spatial symmetry of the MOs
METHOD=none              # break_spatial_sym, fb_loc, pm_loc, none
SET_MO_CLASS=$(qp set_mo_class -d [] -a [] -v [] -i [] -c []) # you can choose the MOs for the localization  
RESET_MO_CLASS=$(qp set_mo_class -c [] -a []) # to reset the core and active MOs after the localization
# you need to know that before launching the script

# Kind of hessian matrix for the orbital optimization
OPT_METHOD=diag          # diag, full
N_DET_MAX_OPT=1e2        # Maximal number of det for the optimization
NB_ITER_OPT=2            # Max number of iteration for the optimization

FILE1=${MOL}_hf
FILE2=${MOL}_no
FILE3=${MOL}_oo
DIR=${MOL}
EZFIO1=${FILE1}.ezfio
EZFIO2=${FILE2}.ezfio
EZFIO3=${FILE3}.ezfio

mkdir $DIR
cd $DIR
cp ../${XYZ}.xyz .

qp_create_ezfio ${XYZ}.xyz -b ${BASIS} -c ${CHARGE} -m ${MULTIPLICITY} -o ${EZFIO1}

echo "### HFOs ###"
echo "Ezfio name:"
echo ${EZFIO1}

qp set_file ${EZFIO1}
qp set determinants read_wf true

# SCF and frozen_core
echo "SCF"
qp run scf > ${FILE1}.scf.out

echo "Frozen core"
qp set_frozen_core > ${FILE1}.fc.out

# Save HFOs
cp -r ${EZFIO1} ${FILE1}_save_mos.ezfio
tar zcf ${FILE1}_save_mos.ezfio.tar.gz ${FILE1}_save_mos.ezfio
rm -r ${FILE1}_save_mos.ezfio

# Selection of the states
if [[ ${SELECT_STATES} == true ]]
then
	echo "State selection 1"
	qp set determinants n_states ${N_MAX_STATES}
	qp set determinants n_det_max ${N_DET_MAX_SELECT}
	qp run fci > ${FILE1}.pre_fci.out
	
	qp edit -s [${STATE_0},${STATE_1}]
fi

# cipsi w HFOs
qp set determinants n_states 2
qp set determinants n_det_max ${N_DET_MAX_HF}
echo "Cipsi with HFOs"
qp run fci > ${FILE1}.fci.out

cp -r ${EZFIO1} ${FILE1}_save_cispi_res.ezfio
tar zcf ${FILE1}_save_cispi_res.tar.gz ${FILE1}_save_cispi_res.ezfio
rm -r ${FILE1}_save_cispi_res.ezfio

# NOs
echo "### NOs ###"
echo "Ezfio name:"
echo ${EZFIO2}

cp -r ${EZFIO1} ${EZFIO2}
qp set_file ${EZFIO2}
echo "Computes NOs"
qp run save_natorb > ${FILE2}.save_natorb.out
qp reset -d

cp -r ${EZFIO2} ${FILE2}_save_mos.ezfio
tar zcf ${FILE2}_save_mos.tar.gz ${FILE2}_save_mos.ezfio
rm -r ${FILE2}_save_mos.ezfio

cp -r ${EZFIO2} ${EZFIO3}

# cipsi w NOs
# Selection of the states
if [[ ${SELECT_STATES} == true ]]
then
        echo "State selection 2"
        qp set determinants n_states ${N_MAX_STATES}
        qp set determinants n_det_max ${N_DET_MAX_SELECT}
        qp reset -d
        qp run fci > ${FILE2}.pre_fci.out

        qp edit -s [${STATE_0},${STATE_1}]
fi

qp set determinants n_states 2
qp set determinants n_det_max ${N_DET_MAX_NO}
echo "Cipsi with NOs"
qp run fci > ${FILE2}.nofci.out

cp -r ${EZFIO2} ${FILE2}_save_cipsi_res.ezfio
tar zcf ${FILE2}_save_cipsi_res.tar.gz ${FILE2}_save_cipsi_res.ezfio
rm -r ${FILE2}_save_cipsi_res.ezfio

echo "### OOs ###"
echo "Ezfio name:"
echo ${EZFIO3}

qp set_file ${EZFIO3}

if [[ ${METHOD} == break_spatial_sym ]] 
then
	echo "Spatial symmetry breaking"
	qp set orbital_optimization security_mo_class false
	qp set orbital_optimization angle_pre_rot 1e-3
	${SET_MO_CLASS}  #qp set_mo_class -d [] -a [] -v []
	qp run break_spatial_sym > ${FILE3}.break_sym.out
    ${RESET_MO_CLASS} #qp set_mo_class -c [] -a []

    cp -r ${EZFIO3} ${FILE3}_save_brk_mos.ezfio
	tar zcf ${FILE3}_save_brk_mos.tar.gz ${FILE3}_save_brk_mos.ezfio
    rm -r ${FILE3}_save_brk_mos.ezfio

elif [[ ${METHOD} == fb_loc ]]
then
	echo "Foster-Boys localization"
    qp set orbital_optimization localization_method boys
    qp set orbital_optimization localization_max_nb_iter 1e4
	qp set orbital_optimization angle_pre_rot 1e-3
    ${SET_MO_CLASS} #qp set_mo_class -d [] -a [] -v []
    qp run localization > ${FILE3}.loc.out
    ${RESET_MO_CLASS} #qp set_mo_class -c [] -a []

    cp -r ${EZFIO3} ${FILE3}_save_fb_mos.ezfio
	tar zcf ${FILE3}_save_fb_lo_mos.tar.gz ${FILE3}_save_fb_mos.ezfio
    rm -r ${FILE3}_save_fb_mos.ezfio

elif [[ ${METHOD} == pm_loc ]]
then
	echo "Pipek-Mezey localization"
	qp set orbital_optimization localization_method pipek
	qp set orbital_optimization localization_max_nb_iter 1e4
	qp set orbital_optimization angle_pre_rot 1e-3
	${SET_MO_CLASS} #qp set_mo_class -d [] -a [] -v []
	qp run localization > ${FILE3}.loc.out
	${RESET_MO_CLASS} #qp set_mo_class -c [] -a []

    cp -r ${EZFIO3} ${FILE3}_save_pm_mos.ezfio
	tar zcf ${FILE3}_save_pm_lo_mos.tar.gz ${FILE3}_save_pm_mos.ezfio
    rm -r ${FILE3}_save_pm_mos.ezfio

else
	echo "No spatial symmetry breaking"
fi

# OOs
## Selection of the states
if [[ ${SELECT_STATES} == true ]]
then
	echo "State selection 3"
	qp set determinants n_states ${N_MAX_STATES}
	qp set determinants n_det_max ${N_DET_MAX_SELECT}
    qp reset -d
	qp run fci > ${FILE3}.pre_fci.out
	
	qp edit -s [${STATE_0},${STATE_1}]
fi

qp set determinants n_states 2

## OO parameters
qp set orbital_optimization optimization_method ${OPT_METHOD}
qp set orbital_optimization normalized_st_av_weight true
qp set orbital_optimization start_from_wf true
qp set orbital_optimization truncate_wf true
qp set orbital_optimization n_det_start 100
qp set orbital_optimization n_det_max_opt ${N_DET_MAX_OPT}
qp set orbital_optimization targeted_accuracy_cipsi 1e-5
qp set orbital_optimization optimization_max_nb_iter ${NB_ITER_OPT}
qp set orbital_optimization thresh_opt_max_elem_grad 1e-4

## Orbital optimization
echo "Orbital optimization"
qp run optimization > ${FILE3}.opt.out

## Save OOs
qp reset -d
tar zcf ${FILE3}_save_mos.tar.gz ${EZFIO3}

## Selection of the states
if [[ ${SELECT_STATES} == true ]]
then
	echo "State selection 4"
	qp set determinants n_states ${N_MAX_STATES}
	qp set determinants n_det_max ${N_DET_MAX_SELECT}
    qp reset -d
	qp run fci > ${FILE3}.pre_opt_fci.out
	
	qp edit -s [${STATE_0},${STATE_1}]
fi

# cispi w OOs
qp set determinants n_states 2
qp set determinants n_det_max ${N_DET_MAX_OO}
echo "Cispi with OOs"
qp run fci > ${FILE3}.opt_fci.out

cp -r ${EZFIO3} ${FILE3}_save_cipsi_res.ezfio
tar zcf ${FILE3}_save_cipsi_res.tar.gz ${FILE3}_save_cipsi_res.ezfio
rm -r ${FILE3}_save_cipsi_res.ezfio

# Data extraction and extrapolation of the fci/excitation energies
echo "Data extraction"
## with NOs
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/extract_E_cipsi.py -f ${FILE2}.nofci.out
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/extrapolation_fci.py -f ${FILE2}.nofci.out.dat > ${FILE2}.extrapolation_fci.dat
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/cipsi_error.py ${FILE2}.nofci.out.dat > ${FILE2}.cipsi_error.dat

## with OOs
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/extract_E_cipsi.py -f ${FILE3}.opt_fci.out
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/extrapolation_fci.py -f ${FILE3}.opt_fci.out.dat > ${FILE3}.extrapolation_fci.dat
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/cipsi_error.py ${FILE3}.opt_fci.out.dat > ${FILE3}.cipsi_error.dat

# Reminder
#MOs
#HFOs: ${FILE1}_save_mos.ezfio.tar.gz
#NOs:  ${FILE2}_save_mos.ezfio.tar.gz
#OOs:  ${FILE3}_save_mos.ezfio.tar.gz
#
#State selection
#State selection for cipsi with HFOs:      ${FILE1}.pre_fci.out
#State selection for cipsi with NOs:       ${FILE2}.pre_fci.out
#State selection for orbital optimization: ${FILE3}.pre_fci.out
#State selection for cipsi with OOs:       ${FILE3}.pre_opt_fci.out
#
#Results of CIPSI calculations
#CIPSI with HFOs: ${FILE1}.fci.out
#CIPSI with NOs:  ${FILE2}.fci.out
#CIPSI with OOs:  ${FILE3}.opt_fci.out
#
#Data from the CIPSI calculations
#CIPSI with HFOs: ${FILE1}.fci.out.dat
#CIPSI with NOs:  ${FILE2}.fci.out.dat
#CIPSI with OOs:  ${FILE3}.opt_fci.out.dat
#
#FCI extrapolation from the CIPSI calculations
#CIPSI with NOs: ${FILE2}.extrapolation_fci.dat
#CIPSI with OOs: ${FILE3}.extrapolation_fci.dat
#
#Extrapolation of the excitation energy
#CIPSI with NOs: ${FILE2}.cispi_error.dat
#CIPSI with OOS: ${FILE3}.cispi_error.dat

