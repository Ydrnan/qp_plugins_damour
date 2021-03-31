Doc

Main program : orb_opt
- methods for the hessian : 
* diagonal hessian 
* full hessian 

- methods for the Newton method :
* Normal
* trust region -> ensure the convergence of the algorithm
* cyrus

Debug :
constants.h contains a global variable 'debug'. 
When debug = .True. the program will display lot of informations.
You may modified some things in the different files to see what you want.
For example, create a new variable Iwant = .True. in constants.h, put debug = .False.

Optimization process :
1. Compute the gradient
     From Anderson et. al. (2014) 
     The Journal of Chemical Physics 141, 244104 (2014); doi: 10.1063/1.4904384

2. Compute the hessian matrix
     - Full hessian
         * From Anderson et. al. (2014)
           The Journal of Chemical Physics 141, 244104 (2014); doi: 10.1063/1.4904384
     - Diagonal hessian
         * Same as full hessian but in a way that only compute the diagonal terms
           (and the terms link to diagonal terms by a permutation
           p <-> q, r <-> s, p <-> q  r <-> s)

3. Compute the inverse of the hessian matrix : Hm1
     - In the case of diagonal hessian, Hm1 = diag(1/diagonal elements)
     - In the case of full hessian :
       * Diagonalization of H
       * Compute diag(1/eigenvectors)
       * Back transform to the initial basis
     - Nothing in the case of trust region

4. Compute the product Hm1.g (Hm1g) and transform Hm1.g in a antisymmetric matrix, m_Hm1g
     - Newton method
       * Compute Hm1.g -> vector of size mo_num*(mo_num-1)/2 corresponding to 
         lower/upper diagonal matrix elements
       * Transform the mo_num*(mo_num-1)/2 vector in a antisymmetric matrix, m_Hm1g
     - With trust region
       * Compute the ratio rho for the model quality
       * Modification of the trust region radius in function of rho
       * Compute the step for the next iteration
         > In function of rho, compute the next step or cancel the previous step
         > If the norm of Hm1g <= trust radius, unconstraint solution, lambda = 0.0
         > If the norm of Hm1g > trust radius, Newton method to find the right 
           Lagrange multiplier to put the constraint on the step size

           The Lagrange multiplier is find by solving :
           ||p(lambda)||^2 - trust_radius = 0
           With p(lambda) = Hm1g(lambda)     
         > if the step is rejected : 
           read the previous step and the actual step becomes -previous step                                           

     - With Umrigar method
       * Compute a factor f_t 
       * Compute f_t . Hm1g

5. Compute the rotation matrix R from the matrix m_Hm1g
   - R = exp(m_Hm1g)

6. Compute the new MOs from the previous one :
     New MOs = Old MOs . R



orb_opt : main program

gradient : optimized gradient subroutine
first_gradient : initial gradient subroutine
orb_opt_debug_grad : program to debug the gradient, compute the difference bewteen gradient and first_gradient

hess : optimized full hessian subroutine
first_hess : initial full hessian subroutine
diag_hess : optimized diagonal hessian subroutine
first_diag_hess : initial diagonal hessian subroutine
orb_opt_debug : program to debug the hessian et diagonal hessian, compute the difference between the first and the optimized version

dm_inversion :  subroutine to compute the inverse of the hessian matrix
dm_Hm1g : subroutine to compute the product Hm1.g and the matrix m_Hm1g

trust_region : subroutine for the trust region
trust_newton : subroutine to find the optimal Lagrange multiplier withe the Newton method
dn_e_model : subroutine to compute the predicted energy
dn_rho_model : subroutine to compute rho, the quality of the model
init_nb_iteration : program to initialize the iteration number

test_cyrus : subroutine to compute the factor f_t 
init_cyrus : program to initialize the cyrus method

dm_rotation : subroutine to compute a rotation matrix

dm_newton_test : subroutine to compute the new MOs

in_mat_vec_index : subroutine to compute the index of a vector element from matrix indexes
in_vec_mat : subroutine to compute the indexes of a matrix element from a vector index

dv_mat_to_vec : subroutine to transform a antisymmetric matrix mo_num by mo_num in a vector of size mo_num.(mo-num-1)/2
dm_vec_to_mat : subroutine to transform a vector of size mo_num.(mo-num-1)/2 in a lower/upper diagonal matrix mo_num by mo_num
The in_mat_vec_index and in_vec_mat is equivalent to dv_mat_to_vec and dm_vec_to_mat, the two first ones compute the indexes of the elements

Calculation procedure :

For the moment the orbital optimization is done using a bash script.
The algorithm orb_opt just produce one step of the Newton method.

The "with hands procedure" is the following :
- Create an ezfio directory
- Do a HF calculation
- Frozen core ?
- Set the number of determinants
- Do a CIPSI calculation
- (Run the init_nb_iteration for trust region and cyrus method)
- (Run the init_cyrus for cyrus method)
- Run the orb_opt algorithm
- Diagonalize the Hamiltonian
- Run the orb_opt algorithm
- Diagonalize the Hamiltonian
- ...

This is the simplest procedure for a fix number of determinants.

The calculation will destroy the starting orbitals, so it is recommended to copy the ezfio directory
containing your orbitals.

At the beginning of the orb_opt output, you can find 3 lines with the details of the method :
diagonal hessian, full hessian, trust region, Umrigar method, ...

Example of a script for the previous procedure :

Before the calculation :

qp create_ezfio h2co_cc_pvdz.xyz -b cc-pvdz
qp run scf > h2co_cc_pvdz.ezfio/h2co_cc_pvdz.scf
qp set_frozen_core > h2co_cc_pvdz.ezfio/h2co_cc_pvdz.frozen
qp set determinants n_det_max X 
qp run fci > h2co_cc_pvdz.ezfio/h2co_cc_pvdz.fci

Script :

#!/bin/bash

source /home/ydamour/qp2/quantum_package.rc #!!! The path must be changed !!!
module load intel/2019.0
module load python/3.7.6-gcc-9.2.0
module load gcc/8.2.0

XYZ=h2co
BASIS=cc_pvdz_hf_h_diag_trust
DIR=${XYZ}_${BASIS}.ezfio
FILE=${XYZ}_${BASIS}

# I work in two different directories :
# one for the SCF/CIPSI calculation
# one for the optimiation
# It's a dangerous way because of the path between the directories

PATH_CIPSI=../../../y_calculs # !!! The path must be changed !!!
PATH_OPT=../plugins/qp_plugins_damour/test_newton # !!! The path must be changed !!!

cp -r h2co_cc_pvdz.ezfio ${DIR}

qp set_file $DIR
qp reset -d
qp set determinants read_wf true
qp set determinants mo_label MCSCF
qp set mo_basis mo_label MCSCF

qp set determinants n_det_max XXX
qp run fci > h2co_cc_pvdz.ezfio/h2co_cc_pvdz_XXX.fci
cd ${PATH_OPT}
qp_run init_nb_iteration ${PATH_CIPSI}/${DIR}
#qp_run init_cyrus ${PATH_CIPSI}/${DIR}
cd ${PATH_CIPSI}

for ((i=1 ; 101 - $i ; i++))
do
cd ${PATH_OPT}
qp_run orb_opt ${PATH_CIPSI}/${DIR} > ${PATH_CIPSI}/${DIR}/opt_trash_XXX_${i}.dat
cd ${PATH_CIPSI}
qp run diagonalize_h > ${DIR}/diag_XXX_${i}.dat
grep "N_det" diag_XXX_${i}.dat >> ${DIR}/${FILE}_n_det.dat
grep "Energy of state" diag_XXX_${i}.dat >> ${DIR}/${FILE}_energy.dat
done

paste ${DIR}/${FILE}_n_det.dat ${DIR}/${FILE}_energy.dat > ${DIR}/${FILE}_result.dat

The same script in one directory is :

#!/bin/bash

source /home/ydamour/qp2/quantum_package.rc #!!! The path must be changed !!!
module load intel/2019.0
module load python/3.7.6-gcc-9.2.0
module load gcc/8.2.0

XYZ=h2co
BASIS=cc_pvdz_hf_h_diag_trust
DIR=${XYZ}_${BASIS}.ezfio
FILE=${XYZ}_${BASIS}

cp -r h2co_cc_pvdz.ezfio ${DIR}

qp set_file $DIR
qp reset -d
qp set determinants read_wf true
qp set determinants mo_label MCSCF
qp set mo_basis mo_label MCSCF

qp set determinants n_det_max XXX
qp run fci > h2co_cc_pvdz.ezfio/h2co_cc_pvdz_XXX.fci
qp run init_nb_iteration
#qp run init_cyrus

for ((i=1 ; 101 - $i ; i++))
do
cd ${PATH_OPT}
qp run orb_opt > ${DIR}/opt_trash_XXX_${i}.dat
cd ${PATH_CIPSI}
qp run diagonalize_h > ${DIR}/diag_XXX_${i}.dat
grep "N_det" diag_XXX_${i}.dat >> ${DIR}/${FILE}_n_det.dat
grep "Energy of state" diag_XXX_${i}.dat >> ${DIR}/${FILE}_energy.dat
done

paste ${DIR}/${FILE}_n_det.dat ${DIR}/${FILE}_energy.dat > ${DIR}/${FILE}_result.dat

Optimize the orbitals for one number of determinants is not the most effective way.
The idea is to increase the number of determinants to optimize the orbitals for
different number of determinants.

#!/bin/bash

source /home/ydamour/qp2/quantum_package.rc #!!! The path must be changed !!!
module load intel/2019.0
module load python/3.7.6-gcc-9.2.0
module load gcc/8.2.0

XYZ=h2co
BASIS=cc_pvdz_hf_h_diag_trust
DIR=${XYZ}_${BASIS}.ezfio
FILE=${XYZ}_${BASIS}

cp -r h2co_cc_pvdz.ezfio ${DIR}

qp set_file $DIR
qp reset -d
qp set determinants read_wf true
qp set determinants mo_label MCSCF
qp set mo_basis mo_label MCSCF

Ndet=10

while (Ndet < 100000) do

	qp set determinants n_det_max Ndet
	qp run fci > ${DIR}/${FILE}_${Ndet}.fci

	grep "Summary at N_det = " ${DIR}/${FILE}_${Ndet}.fci >> ${DIR}/${FILE}_n_det.dat
	grep "# E " ${DIR}/${FILE}_${Ndet}.fci >> ${DIR}/${FILE}_energy.dat
	grep "# PT2 " ${DIR}/${FILE}_${Ndet}.fci >> ${DIR}/${FILE}_pt2.dat
	grep "# rPT2 " ${DIR}/${FILE}_${Ndet}.fci >> ${DIR}/${FILE}_rpt2.dat

	qp run init_nb_iteration
	#qp run init_cyrus
	
	for ((i=1 ; 101 - $i ; i++))
	do
		cd ${PATH_OPT}
		qp run orb_opt > ${DIR}/opt_trash_${Ndet}_${i}.dat
		cd ${PATH_CIPSI}
		qp run diagonalize_h > ${DIR}/diag_${Ndet}_${i}.dat
		grep "N_det" diag_${N_det}_${i}.dat >> ${DIR}/${FILE}_n_det_opt.dat
		grep "Energy of state" diag_${Ndet}_${i}.dat >> ${DIR}/${FILE}_energy_opt.dat
		grep "N_det" diag_${N_det}_${i}.dat > ${DIR}/tmp_n_det
	done
	
	qp run pt2 > ${DIR}/${FILE}_${N_det}.pt2 

	grep "Summary at N_det = " ${DIR}/${FILE}_${Ndet}.pt2 >> ${DIR}/${FILE}_n_det.dat
	grep "# E " ${DIR}/${FILE}_${Ndet}.pt2 >> ${DIR}/${FILE}_energy.dat
	grep "# PT2 " ${DIR}/${FILE}_${Ndet}.pt2 >> ${DIR}/${FILE}_pt2.dat
	grep "# rPT2 " ${DIR}/${FILE}_${Ndet}.pt2 >> ${DIR}/${FILE}_rpt2.dat
	
	cat ${DIR}/tmp_n_det | while read Ndet ; do
		Ndet=$(($Ndet*2))
		qp set determinants n_det_max Ndet
	done

done

qp set determinants n_det_max 1e6
qp run fci > ${DIR}/${FILE}_1e6.fci

grep "Summary at N_det = " ${DIR}/${FILE}_1e6.fci >> ${DIR}/${FILE}_n_det.dat
grep "# E " ${DIR}/${FILE}_1e6.fci >> ${DIR}/${FILE}_energy.dat
grep "# PT2 " ${DIR}/${FILE}_1e6.fci >> ${DIR}/${FILE}_pt2.dat
grep "# rPT2 " ${DIR}/${FILE}_1e6.fci >> ${DIR}/${FILE}_rpt2.dat

paste ${DIR}/${FILE}_n_det_opt.dat ${DIR}/${FILE}_energy_opt.dat > ${DIR}/${FILE}_result_opt.dat

paste ${DIR}/${FILE}_n_det.dat ${DIR}/${FILE}_energy.dat > ${DIR}/${FILE}_tmp1.dat
paste ${DIR}/${FILE}_tmp1.dat ${DIR}/${FILE}_pt2.dat > ${DIR}/${FILE}_tmp2.dat
paste ${DIR}/${FILE}_tmp2.dat ${DIR}/${FILE}_rpt2.dat > ${DIR}/${FILE}_result.dat

For a large number of determinants the optimization of orbitals has a very weak impact on the energy.
The biggest part of the orbital optimization is done for the small number of determinants. 
Clearly it depends of the molecule and the basis.
Since there is no convergence criterion the number of step for each optimization must be changed.

!!! Warning !!!

For the moment some informations are store in temporary files...
So for trust region or cyrus method it is not possible to run more than one calculation 
at the same time...

Trust region is mandatory to optimize orbitals. 
Trust region ensure the convergence of the optimization.





