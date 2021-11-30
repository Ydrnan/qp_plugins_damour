Installation on QP2
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

To optimize the orbitals for this number of determinants:
qp set determinants read_wf true
qp set determinants mo_label MCSCF
qp set mo_basis mo_label MCSCF
qp set ao_two_e_erf_ints io_ao_two_e_integrals_erf write
qp run org_orb_opt_trust_v2 > optimization.dat

The ezfio directory now contains the optimized MOs

Some parameters can be changed with qp edit in the Orbital_optimization section

If you modify the .org files, don't forget to do
./TANGLE_org_mode.sh
ninja

The documentation can be read using:
Ctrl-C Ctrl-e l p
after opening the filename.org in emacs. It will produce a 
filename.pdf. 
!!! Warning: the documentation can contain some errors !!!

============================================
Further improvements:
- Cleaner repo
- Correction of the errors in the documentations
- Elliptical trust region
- Quasi Newton
