# Selected CI
Cipsi with some constraints on the excitation degree or seniority.  
  
Note:  
All these kinds of calculation can be done by setting by hands  
the desired maximal excitation / seniority such as:
```
qp set cipsi excitation_max #a number here
qp set cipsi seniority_max #a number here
```

## Excitation-based CI, CIPSI version
CISD, CISDT, CISTQ  
```
qp run cisd_cipsi
qp run cisdt_cipsi
qp run cisdtq_cipsi
```

## Excitation-based CI, non CIPSI version
Only for the same number of electron alpha and beta for the moment.   
CISDT, CISTQ  
```
qp run cisdt
qp run cisdtq
```

## Seniority-based CI, CIPSI version
sCI0, sCI2, sCI4, sCI6  
```
qp run s_ci0_cipsi
qp run s_ci2_cipsi
qp run s_ci4_cipsi
qp run s_ci6_cipsi
```

## Hierarchy based CI, CISPI version
hCI-1, hCI-1.5, hCI-2, hCI-2.5, hCI-3  
```
qp run hci_1_cipsi.irp.f
qp run hci_1d5_cipsi.irp.f
qp run hci_2_cipsi.irp.f
qp run hci_2d5_cipsi.irp.f
qp run hci_3_cipsi.irp.f
```

## PT2
The pt2 calculations done during these calculations are only done in  
the considered Hilbert subspace. But a standart pt2 calculation can be  
done with:
```
qp run pt2
```

# Settings
Don't forget to set:  
- the maximal number of determinants:  
```
qp set determinants n_det_max #a number here
```
- the number of states:  
```
qp set determinants n_states #a number here
```
  
And possibly:  
- if you run only on one node:  
```
qp set davidson distributed_davidson false
```
- to speed up/ reduced the memory:  
```
qp set davidson csf_based true
```
- if you want to avoid pt2 calculation:  
```
qp set perturbation do_pt2 false
```
- the cipsi selection stops after reaching a given PT2,  
you can change this with
```
qp set perturbation pt2_max #a number here
```

