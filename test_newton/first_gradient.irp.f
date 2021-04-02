subroutine first_gradient(n,v_grad)

  include 'constants.h'

  implicit none

  !===================================================================
  ! Compute the gradient of energy with respects to orbital rotations
  !===================================================================

  ! Check if read_wf = true, else :
  ! qp set determinant read_wf true

  END_DOC

  ! in
  integer, intent(in) :: n
  ! n : integer, n = mo_num*(mo_num-1)/2
  
  ! out
  double precision, intent(out) :: v_grad(n)
  ! v_grad : double precision vector of length n containeing the gradient

  ! internal
  double precision, allocatable :: grad(:,:),A(:,:)
  double precision :: norm
  integer :: i,p,q,r,s,t
  integer :: istate
  ! grad : double precision matrix containing the gradient before the permutation
  ! A : double precision matrix containing the gradient after the permutation
  ! norm : double precision number, the norm of the vector gradient
  ! i,p,q,r,s,t : integer, indexes 
  ! istate : integer, the electronic state

  ! Function
  double precision :: get_two_e_integral, norm2
  ! get_two_e_integral :  double precision function that gives the two e integrals
  ! norm2 : double precision function that gives the norm of a vector
 
  ! Provided :
  ! mo_one_e_integrals : mono e- integrals
  ! get_two_e_integral : two e- integrals
  ! one_e_dm_mo_alpha, one_e_dm_mo_beta : one body density matrix
  ! two_e_dm_mo : two body density matrix

  !============
  ! Allocation
  !============

  allocate(grad(mo_num,mo_num),A(mo_num,mo_num))

  !=============
  ! Calculation
  !=============

  if (debug) then
    print*,'Enter in first_gradient'
  endif

  v_grad = 0d0

  !do istate = 1, N_states
  istate = 1
    do p = 1, mo_num
      do q = 1, mo_num
         grad(p,q) = 0d0
         do r = 1, mo_num
           grad(p,q) = grad(p,q) + mo_one_e_integrals(p,r) &
                          * (one_e_dm_mo_alpha(r,q,istate) + one_e_dm_mo_beta(r,q,istate)) &
                         - mo_one_e_integrals(r,q) &
                          * (one_e_dm_mo_alpha(p,r,istate) + one_e_dm_mo_beta(p,r,istate))

        enddo

       do r = 1, mo_num
         do s = 1, mo_num
           do t= 1, mo_num

           grad(p,q) = grad(p,q) &
                   + get_two_e_integral(p,t,r,s,mo_integrals_map) * two_e_dm_mo(r,s,q,t) &
                   - get_two_e_integral(r,s,q,t,mo_integrals_map) * two_e_dm_mo(p,t,r,s)
           enddo
          enddo
        enddo
      enddo
    enddo
  !enddo

 ! Debug ocaml
!  print*,'two e rdm'
!  do p=1,mo_num
!    do q=1,mo_num
!      do r=1,mo_num
!        do s=1,mo_num
!          print*,p,q,r,s,two_e_dm_mo(p,q,r,s,1)
!        enddo
!      enddo
!    enddo
!  enddo
!
!  print*,'bi int'
!  do p=1,mo_num
!    do q=1,mo_num
!      do r=1,mo_num
!        do s=1,mo_num
!          if (ABS(get_two_e_integral(p,q,r,s,mo_integrals_map))>1.d-14) then
!            print*,p,q,r,s, get_two_e_integral(p,q,r,s,mo_integrals_map)
!          else
!            print*,p,q,r,s,0d0
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!
!  print*,'mono int'
!  do p=1,mo_num
!    print*,mo_one_e_integrals(p,:)
!  enddo
!
!  print*, 'one e rdm'
!  do p=1,mo_num
!    print*, (one_e_dm_mo_alpha(p,:,istate) + one_e_dm_mo_beta(p,:,istate))
!  enddo

 !Ecriture des intégrales dans un fichier pour le lire avec OCaml

  if (ocaml) then
    open(unit=10,file='../../../../../../App_y/miniconda3/Work_yann/one_e_dm.dat')
    do p = 1, mo_num
      do q = 1, mo_num
        write(10,*) p, q, (one_e_dm_mo_alpha(p,q,istate) + one_e_dm_mo_beta(p,q,istate))
      enddo
    enddo
    close(10)
  
    open(unit=11,file='../../../../../../App_y/miniconda3/Work_yann/one_e_int.dat')
    do p = 1, mo_num
      do q = 1, mo_num
        write(11,*) p, q, mo_one_e_integrals(p,q)
      enddo
    enddo
    close(11)
    
    open(unit=12,file='../../../../../../App_y/miniconda3/Work_yann/two_e_int.dat')
    do p = 1, mo_num
      do q = 1, mo_num
        do r = 1, mo_num
          do s = 1, mo_num
            write(12,*) p, q, r, s, get_two_e_integral(p,q,r,s,mo_integrals_map)
          enddo
        enddo
      enddo
    enddo
    close(12)
    
    open(unit=13,file='../../../../../../App_y/miniconda3/Work_yann/two_e_dm.dat')
    do p = 1, mo_num
      do q = 1, mo_num
        do r = 1, mo_num
          do s = 1, mo_num
            write(13,*) p, q, r, s, two_e_dm_mo(p,q,r,s)
          enddo
        enddo
      enddo
    enddo
    close(13)
  endif

 ! Debug ocaml
!  double precision, allocatable :: two_e_integrals(:,:,:,:)
!  allocate(two_e_integrals(mo_num,mo_num,mo_num,mo_num))
!
!  do p=1,mo_num
!    do q=1,mo_num
!      do r=1,mo_num
!        do s=1,mo_num
!          two_e_integrals(p,q,r,s) = get_two_e_integral(p,q,r,s,mo_integrals_map)
!        enddo
!      enddo
!    enddo
!  enddo
!
!
!  print*,'two e int'
!  do p = 1,mo_num
!    do q= 1, mo_num
!      write(*,'(100(F10.5))') two_e_integrals(p,q,:,:)
!    enddo
!  enddo
!
!  print*,'two e rdm'
!  do p = 1,mo_num
!    do q= 1, mo_num
!      write(*,'(100(F10.5))') two_e_dm_mo(p,q,:,:,1)
!    enddo
!  enddo

  ! Conversion mo_num*mo_num matrix to mo_num(mo_num-1)/2 vector
  do i=1,n
    call in_mat_vec_index(i,p,q)
    v_grad(i)=(grad(p,q) - grad(q,p))
  enddo  

  ! Display, vector containing the gradient elements 
  if (debug) then  
    print*,'Vector containing the gradient :'
    write(*,'(100(F10.5))') v_grad(1:n)
  endif  

  ! Norm of the vector
  norm = norm2(v_grad)
  print*, 'Norm : ', norm

  ! Matrix gradient
  A = 0d0
  do q=1,mo_num
    do p=1,mo_num
      A(p,q) = grad(p,q) - grad(q,p)
    enddo
  enddo

  ! Display, matrix containting the gradient elements
  if (debug) then
    print*,'Matrix containing the gradient :'
    do i = 1, mo_num
      write(*,'(100(F10.5))') A(i,1:mo_num)
    enddo
  endif

  !==============
  ! Deallocation
  !==============

  deallocate(grad,A)

  if (debug) then
    print*,'Leaves first_gradient'
  endif

end subroutine
