# Plugin to compute dipole moment and oscillator strength

## CIS
Like the original cis:  
```
qp run cis_w_dipoles
```

## CISD
Like the original cisd:  
```
qp run cisd_w_dipoles
```

## FCI
Like the original fci:  
```
qp run fci_w_dipoles
```

## Results
- Dipole in au and Debye  
- Transition dipole between the ground state and the excited states  
  in au and Debye.  
- Oscillator strength between the ground state and the excited states  
  in length gauge.  

If you want to print the transition dipoles and the oscillator strengths   
for all the transitions, set:
```
qp set dipoles print_all_transitions true
```

If you want to print the n_det_print first determinants of each state:  
```
qp set dipoles print_det_state true
```

