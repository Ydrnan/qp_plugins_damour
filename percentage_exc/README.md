==============
percentage_exc
==============

# Installation
```
cd $QP_ROOT/plugins
qp plugins install percent_exc
cd qp_plugins_damour/percent_exc
./TANGLE_org_mode.sh
ninja
```

# Extraction of the CI coefficients from CC amplitudes  
Run a cc calculation and write the T1 and T2 amplitudes disk
```
qp set cc_utils cc_write_t1 true   
qp set cc_utils cc_write_t2 true
qp run ccsd_spin_orb
```
Then the CI coefficients can be extracted
```
qp run extract_c
```
The C3 and C4 coefficients can be also extracted as product of t1 and t2  
by setting
```
qp set percentage_exc extract_c3 true  
qp set percentage_exc extract_c4 true
```
But the their number can rapidly blow up. In addition, the number  
of CI coefficients can be limited by playing on the threshold   
of the selection
```
qp set percentage_exc thresh_extract_c 1e-4 # !!!In intermediate normalization!!!
```

# Percentage of excitations
Format exponential:  
```
qp set percent_exc percentage_in_exp true
```
Percentage C:   
```
qp run print_percent_c
```
Percentage T:
```
qp run print_percent_t
```

## Number of excitation
```
qp set percent_exc nb_excitation_max
```

For the percentage C, a non-HF reference can be chosen:  
```
qp set percent_exc ref_type 2 ! 1 for HF ref
```
To do that you have to create a file "ref.txt" containing on the first  
line the number of determinants for your reference wave function, and one  
determinant per line by writing the 2 N_int integer of each  
determinants (the first N_int integer for the alpha part and the last  
N_int integer for the beta part).
