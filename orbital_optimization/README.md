Installation on qp2   
Go to the qp2 directory 
``` 
./bin/qpsh  
cd plugins  
git clone -b dev https://github.com/Ydrnan/qp_plugins_damour  
qp_plugins install orbital_optimization  
cd qp_plugins_damour/trust_region
./TANGLE_org_mode.sh  
ninja
cd ../orbital_optimization  
./TANGLE_org_mode.sh  
ninja  
``` 
Please, use the ifort compiler  
  
Some parameters can be changed with qp edit in the Orbital_optimization section 
 
If you modify the .org files, don't forget to do:  
``` 
./TANGLE_org_mode.sh  
ninja  
```  

The documentation can be read using:  
Ctrl-C Ctrl-e l p  
after opening the filename.org in emacs. It will produce a  
filename.pdf.  
(Not available for all the files)  
!!! Warning: the documentation can contain some errors !!! 


Please, before the orbital optimization run:  
``` 
qp set determinants read_wf true  
qp set ao_two_e_erf_ints io_ao_two_e_integrals_erf write  
``` 

Different methods are available:  
- full hessian  
``` 
qp set orbital_optimization optimization_method full  
```  
- diagonal hessian  
``` 
qp set orbital_optimization optimization_method diag  
``` 
- identity matrix  
``` 
qp set orbital_optimization optimization_method none  
``` 

After the optimization the ezfio contains the optimized orbitals
 
## For a fixed number of determinants
To optimize the MOs for the actual determinants:  
``` 
qp run orb_opt
``` 
 
## For a complete optimization, i.e, with a larger and larger wave function
To optimize the MOs with a larger and larger wave function:  
``` 
qp run optimization  
``` 

The results are stored in the file "result_opt.dat",
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
  
### End of the optimization
You can choos the number of determinants after what the 
optimization will stop:
```
qp set orbital_optimization n_det_max_opt 1e5 # or any number > n_det_start 
```
## Weight of the states
You can change the weights of the differents states directly in qp edit.  
It will affect ths weights used in the orbital optimization.

# Tests
To run the tests:  
``` 
qp test
``` 
