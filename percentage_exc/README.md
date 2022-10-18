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
