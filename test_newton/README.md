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
dn_e_model : subroutine to compute the predicted energy
dn_rho_model : subroutine to compute rho, the quality of the model

rotation_matrix : subroutine to compute a rotation matrix
rotation_matrix_omp : subroutine to compute a rotation matrix (omp version)

apply_mo_rotation : subroutine to compute the new MOs
apply_mo_rotation_omp : subroutine to compute the new MOs (omp version)

in_mat_vec_index : subroutine to compute the index of a vector element from matrix indexes
in_vec_mat : subroutine to compute the indexes of a matrix element from a vector index

dv_mat_to_vec : subroutine to transform a antisymmetric matrix mo_num by mo_num in a vector of size mo_num.(mo-num-1)/2
dm_vec_to_mat : subroutine to transform a vector of size mo_num.(mo-num-1)/2 in a lower/upper diagonal matrix mo_num by mo_num
The in_mat_vec_index and in_vec_mat_index is equivalent to dv_mat_to_vec and dm_vec_to_mat, the two first ones compute the indexes of the elements

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

Example of a script for the previous procedure :

For a large number of determinants the optimization of orbitals has a very weak impact on the energy.
The biggest part of the orbital optimization is done for the small number of determinants. 
Clearly it depends of the molecule and the basis.

!!! Warning !!!

Trust region is mandatory to optimize orbitals. 
Trust region ensure the convergence of the optimization.





