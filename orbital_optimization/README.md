iInstallation on QP2 
Go to the qp2 directory 
 
./bin/qpsh  
cd plugins  
git clone -b dev https://github.com/Ydrnan/qp_plugins_damour  
qp_plugins install orbital_optimization  
cd qp_plugins_damour/orbital_optimization  
./TANGLE_org_mode.sh  
cd $QP_ROOT  
ninja  
  
Please, use the ifort compiler  
  
Some parameters can be changed with qp edit in the Orbital_optimization section 
 
If you modify the .org files, don't forget to do:  
./TANGLE_org_mode.sh  
ninja  
 
The documentation can be read using: 
Ctrl-C Ctrl-e l p 
after opening the filename.org in emacs. It will produce a  
filename.pdf.  
(Not available for all the files) 
!!! Warning: the documentation can contain some errors !!! 

# Orbital localisation
To localize the MOs: 
qp run localization 
 
After that the ezfio directory contains the localized MOs 
 
But the mo_class must be defined before, run qp set_mo_class -q for more information 
 
## Foster-Boys & Pipek-Mezey
Foster-Boys:  
qp set orbital_optimization localization_method boys 
 
Pipek-Mezey:  
qp set orbital_optimization localization_method pipek 

# Orbital optimization 
Please before the orbital optimization run:  
qp set determinants read_wf true  
qp set ao_two_e_erf_ints io_ao_two_e_integrals_erf write  

Different methods are available:  
- full hessian  
qp set orbital_optimization optimization_method full  
 
- diagonal hessian (not available with org_orb_opt_trust_v2)  
qp set orbital_optimization optimization_method diag  
 
- identity matrix  
qp set orbital_optimization optimization_method none  
 
## For a fixed number of determinants
To optimize the MOs for the actual determinants:  
qp run org_orb_opt_trust  
or  
qp run org_orb_opt_trust_v2 (to compute the hessian only for the active MOs)  
 
## For a complete optimization, i.e, with a larger and larger wave function
To optimize the MOs with a larger and larger wave function:  
qp run optimization  
 
The results are stored in the file "result_orbital_optimization.dat",
with the following format:  
(1) (2) (3) (4)  
1: Number of determinants in the wf,  
2: Cispi energy before the optimization,   
3: Cipsi energy after the optimization,  
4: Energy difference between (2) and (3).  
 
The optimization process if the following: 
- we do a first cipsi step to obtain a small number of determinants in the wf 
- we run an orbital optimization for this wf 
- we do a new cipsi step to double the number of determinants in the wf 
- we run an orbital optimization for this wf 
- ... 
- we do that until the energy difference between (2) and (3) is  
  smaller than the targeted accuracy for the cispi (targeted_accuracy_cipsi in qp edit) 
  or the wf is larger than a given size (n_det_max_opt in qp_edit) 
- after that you can reset your determinants (qp reset -d) and run a clean Cispi calculation 
 
============================================  
Further improvements: 
- Cleaner repo 
- Correction of the errors in the documentations 
- Elliptical trust region 
- Quasi Newton 
 
To run the tests:  
./auto_test.sh  
or   
./auto_check.sh  
 
But it will take a long time, so it's better to stop the tests after the tests on the benzene 
===========================================  
