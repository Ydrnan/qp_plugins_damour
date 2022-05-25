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
qp run orb_opt_trust  # old version
``` 
or  
``` 
qp run orb_opt_trust_v2  # new version
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
  
Some settings are available for that.  
### Starting number of determinants
You can choose the number of determinants to start the  
orbital optimization with:
```
qp set orbital_optimization n_det_start 100 # or n_det_max_opt > any number > 0
```
in order to do a first cipsi until this number before  
the optimization.  
You can also choose to start from a given cipsi wave function:
```
qp set orbital_optimization start_from_wf true
```
or to truncate the actual cispi wave function after 
N determinants:
```
qp edit -n N 
```

### End of the optimization
You can choos the number of determinants after what the 
optimization will stop:
```
qp set orbital_optimization n_det_max_opt 1e5 # or any number > n_det_start 
```
You can set the targeted accuracy to stop the optimization when the   
gain in energy for any optimization is smaller than a threshold:
```
qp set orbital optimization targeted_accuracy_cipsi 1e-4
```
Generally, it's better to put the researched accuracy divided by 10,   
because some steps can provide a small gain before a bigger one.

## Weight of the states
You can enforce the weights of the states to be equal and normalized with:
```
qp set orbital_optimization normalized_st_av_weight true
```
or just enforcing the normalization of the weigths:
```
qp set orbital_optimization normalized_weight true # Not tested

```

# Tests
To run the tests:  
``` 
./auto_test_v2.sh  
``` 
 
But it will take a long time, so it's better to stop the tests after  
the tests on the gradient (the tests on the hessians are definitely too long).  
Note: the error on the gradient and hessian must be smaller than 1e-12.  

# Further improvements: 
- Cleaner repo 
- Correction of the errors in the documentations 
- Elliptical trust region 
- Quasi Newton 
 
