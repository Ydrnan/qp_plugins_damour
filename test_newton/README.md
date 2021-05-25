Doc

Main program : orb_opt_trust
- methods for the hessian : 
    * diagonal hessian 
    * full hessian 

Program to extract data :
- extrapolation : program to extract data in the cipsi output and compute the extrapolated correlation energy (mH)
This program creates some files. Paid attention you could encounter format problems with this program.

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

Subroutines :

gradient : optimized gradient subroutine
first_gradient : initial gradient subroutine

hess : optimized full hessian subroutine
first_hess : initial full hessian subroutine
diag_hess : optimized diagonal hessian subroutine
first_diag_hess : initial diagonal hessian subroutine

debug_hessian : program to debug the hessian and diagonal hessian, it computes the difference between the first and the optimized version
debug_gradient : program to debug the gradient, it computes the difference between the first and the optimized version
debug_rotation : program to debug the rotation matrix, it computes the difference between the omp version and the normal version
debug_new_mos : program to debug the new_mos, it computes the difference between the omp version and the normal version

matrix_inversion :  subroutine to compute the inverse of the hessian matrix
Hm1g : subroutine to compute the product Hm1.g and the matrix m_Hm1g

trust_region : subroutine for the trust region
trust_newton : subroutine to find the optimal Lagrange multiplier withe the Newton method
trust_e_model : subroutine to compute the predicted energy
rho_model : subroutine to compute rho, the quality of the model

rotation_matrix : subroutine to compute a rotation matrix
rotation_matrix_omp : subroutine to compute a rotation matrix (omp version)

apply_mo_rotation : subroutine to compute the new MOs
apply_mo_rotation_omp : subroutine to compute the new MOs (omp version)

mat_to_vec_index : subroutine to compute the index of a vector element from a matrix index
vec_to_mat_index : subroutine to compute the index of a matrix element from a vector index

/error/mat_to_vec : subroutine to transform a antisymmetric matrix mo_num by mo_num in a vector of size mo_num.(mo-num-1)/2
/error/vec_to_mat : subroutine to transform a vector of size mo_num.(mo-num-1)/2 in a lower/upper diagonal matrix mo_num by mo_num

Calculation procedure :

The "with hands procedure" is the following :
- Create an ezfio directory
- Do a HF calculation
- Frozen core ?
- Set the number of determinants
- Do a CIPSI calculation
- Run the orb_opt_trust algorithm

This is the simplest procedure for a fix number of determinants.

The calculation will destroy the starting orbitals, so it is recommended to copy the ezfio directory
containing your orbitals.

!!! Warning !!!

Trust region is mandatory to optimize orbitals.
Trust region ensure the convergence of the optimization.

Example of a script for the previous procedure :

For a large number of determinants the optimization of orbitals has a very weak impact on the energy.
The biggest part of the orbital optimization is done for the small number of determinants. 
Clearly it depends of the molecule and the basis.

#!/bin/bash
#SBATCH -p xeonv3 -N 1 -n 1 -c 24 --exclusive
source /home/ydamour/qp2/quantum_package.rc
module load intel/2019.0
module load python/3.7.6-gcc-9.2.0
module load gcc/8.2.0

XYZ=benzene
BASIS=HF_opt_S2
DIR=${XYZ}_${BASIS}.ezfio
FILE=${XYZ}_${BASIS}
INITIAL_DIR

cp -r ${INITIAL_DIR} ${DIR}

qp set_file $DIR
qp reset -a
qp reset -d
qp run scf > ${DIR}/${FILE}.scf

# frozen core ?
qp set_frozen_core > ${DIR}/${FILE}.frozen

qp set determinants read_wf true
qp set determinants mo_label MCSCF
qp set mo_basis mo_label MCSCF
qp set ao_two_e_erf_ints io_ao_two_e_integrals_erf write

# S^2 true or false
#qp set determinants s2_eig true
qp set determinants s2_eig false

# Number of determinants
Ndet=5

while [ ${Ndet} -lt 200000 ]
do
    echo ${Ndet}
    qp set determinants n_det_max ${Ndet}
    qp run fci > ${DIR}/${FILE}_${Ndet}.fci
 
    grep "Summary at N_det = " ${DIR}/${FILE}_${Ndet}.fci | tail -1 >> ${DIR}/${FILE}_n_det.dat
    grep "# E   " ${DIR}/${FILE}_${Ndet}.fci | tail -1 >> ${DIR}/${FILE}_energy.dat
    grep "# PT2   " ${DIR}/${FILE}_${Ndet}.fci | tail -1 >> ${DIR}/${FILE}_pt2.dat
    grep "# rPT2   " ${DIR}/${FILE}_${Ndet}.fci | tail -1 >> ${DIR}/${FILE}_rpt2.dat

    qp run orb_opt_trust > ${DIR}/${FILE}_opt_trash_${Ndet}.dat
    grep "Max element in gardient :" ${DIR}/${FILE}_opt_trash_${Ndet}.dat > ${DIR}/${FILE}_grad_${Ndet}.dat
    grep "Number of negative eigenvalues :"  ${DIR}/${FILE}_opt_trash_${Ndet}.dat > ${DIR}/${FILE}_eval_${Ndet}.dat        
    grep "e_val < 0 :"  ${DIR}/${FILE}_opt_trash_${Ndet}.dat > ${DIR}/${FILE}_eval_value_${Ndet}.dat
    grep "Energy of state    1" ${DIR}/${FILE}_opt_trash_${Ndet}.dat > ${DIR}/${FILE}_opt_energy_${Ndet}.dat

    qp run pt2 > ${DIR}/${FILE}_${Ndet}.pt2 

    grep "Summary at N_det = " ${DIR}/${FILE}_${Ndet}.pt2 | tail -1 >> ${DIR}/${FILE}_n_det.dat
    grep "# E   " ${DIR}/${FILE}_${Ndet}.pt2 | tail -1 >> ${DIR}/${FILE}_energy.dat
    grep "# PT2   " ${DIR}/${FILE}_${Ndet}.pt2 | tail -1 >> ${DIR}/${FILE}_pt2.dat
    grep "# rPT2   " ${DIR}/${FILE}_${Ndet}.pt2 | tail -1 >> ${DIR}/${FILE}_rpt2.dat

Ndet=$[${Ndet}*2]    
done

paste ${DIR}/${FILE}_n_det_opt.dat ${DIR}/${FILE}_energy_opt.dat > ${DIR}/${FILE}_result_opt.dat
paste ${DIR}/${FILE}_n_det.dat ${DIR}/${FILE}_energy.dat > ${DIR}/${FILE}_tmp1.dat
paste ${DIR}/${FILE}_tmp1.dat ${DIR}/${FILE}_pt2.dat > ${DIR}/${FILE}_tmp2.dat
paste ${DIR}/${FILE}_tmp2.dat ${DIR}/${FILE}_rpt2.dat > ${DIR}/${FILE}_result.dat

# CIPSI calculation
DIR2=${FILE}_cipsi.ezfio
FILE2=${FILE}_cipsi
cp -r ${DIR} ${DIR2}

qp set_file ${DIR2}
qp reset -d
qp set determinants read_wf true
qp set determinants mo_label MCSCF
qp set mo_basis mo_label MCSCF

# S^2 true or false
#qp set determinants s2_eig true
qp set determinants s2_eig false

# Number of determinants
qp set determinants n_det_max 3e6

qp run fci > ${DIR2}/${FILE2}.fci

grep "Summary at N_det = " ${DIR2}/${FILE2}.fci >> ${DIR2}/${FILE2}_n_det.dat
grep "# E   " ${DIR2}/${FILE2}.fci >> ${DIR2}/${FILE2}_energy.dat
grep "# PT2   " ${DIR2}/${FILE2}.fci >> ${DIR2}/${FILE2}_pt2.dat
grep "# rPT2   " ${DIR2}/${FILE2}.fci >> ${DIR2}/${FILE2}_rpt2.dat

paste ${DIR2}/${FILE2}_n_det_opt.dat ${DIR2}/${FILE2}_energy_opt.dat > ${DIR2}/${FILE2}_result_opt.dat

paste ${DIR2}/${FILE2}_n_det.dat ${DIR2}/${FILE2}_energy.dat > ${DIR2}/${FILE2}_tmp1.dat
paste ${DIR2}/${FILE2}_tmp1.dat ${DIR2}/${FILE2}_pt2.dat > ${DIR2}/${FILE2}_tmp2.dat
paste ${DIR2}/${FILE2}_tmp2.dat ${DIR2}/${FILE2}_rpt2.dat > ${DIR2}/${FILE2}_result.dat

qp run extrapolation
echo ${DIR2}
echo ${FILE2}
