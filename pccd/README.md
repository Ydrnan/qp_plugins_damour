# pCCD

To run a pCCD calculation:  
```
qp run pccd
```

## Guess
The guess availables are:
- none, i.e., the amplitudes are set to 0
- MP, i.e., the amplitudes are set to the MP2 ones
The guess is chosen with:
```
qp set cc_utils cc_guess_t2 #guess_type
```

## Update
The amplitudes can be updated using the standard approach or the jacobian.  
The diis update is not available yet. The jacobian can be used with:  
```
qp set pccd pccd_update_t2 full
```
for the full jacobian
```
qp set pccd pccd_update_t2 diag
```
for the diagonal jacobian
```
qp set pccd pccd_update_t2 none
```
for the standard update method.
