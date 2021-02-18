subroutine gradient(n,v_grad)
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  
  ! Check if read_wf = true, else :
  ! qp set determinant read_wf true 
 
  END_DOC
  !double precision, allocatable :: grad(:,:)
  !double precision, allocatable :: v_grad(:)
  !double precision, intent(out) :: grad(mo_num,mo_num)

  integer, intent(in) :: n
  double precision, intent(out) :: v_grad(n)
  double precision, allocatable :: grad(:,:)

  double precision :: get_two_e_integral
  double precision :: accu1, accu2
  integer :: i,p,q,r,s,t
  integer :: qorb, porb,rorb,sorb,torb, istate

  !print*, mo_num
  
  !n=mo_num*(mo_num-1)/2

  !============
  ! Allocation
  !============
  !allocate(v_grad(n))
  allocate(grad(mo_num,mo_num))

  !=============
  ! Calculation
  !=============

  v_grad = 0d0
  grad = 0d0

  !do istate = 1, N_states
  istate = 1
    do q = 1, n_act_orb
      qorb = list_act(q)
      do p = 1, n_act_orb
       porb = list_act(p)
       accu1 = 0d0
         do r = 1, n_act_orb
           rorb = list_act(r)

           accu1 = accu1 + mo_one_e_integrals(porb,rorb) &
                          * (one_e_dm_mo_alpha(r,q,istate) + one_e_dm_mo_beta(r,q,istate)) &
                         - mo_one_e_integrals(rorb,qorb) &
                          * (one_e_dm_mo_alpha(p,r,istate) + one_e_dm_mo_beta(p,r,istate)) 
 
!  do q=1,mo_num
!    do p=1,mo_num
!      accu1 = 0d0
!      do r=1,mo_num
!        accu1 = accu1 + mo_one_e_integrals(p,r) &
!                       * (one_e_dm_mo_alpha(r,q,1) + one_e_dm_mo_beta(r,q,1)) &
!                      - mo_one_e_integrals(r,q) &
!                       * (one_e_dm_mo_alpha(p,r,1) + one_e_dm_mo_beta(p,r,1)) !&
                      !- mo_one_e_integrals(q,r) &
                      ! * (one_e_dm_mo_alpha(r,p,1) + one_e_dm_mo_beta(r,p,1)) &
                      !+ mo_one_e_integrals(r,p) &
                      ! * (one_e_dm_mo_alpha(q,r,1) + one_e_dm_mo_beta(q,r,1))
!      enddo
!      grad(p,q) = grad(p,q) +  accu1
!    enddo
!  enddo  
    
        enddo
        grad(p,q) = grad(p,q) + accu1
!      enddo
!    enddo
  !enddo

 
  !do istate = 1, N_states 
!  istate = 1
!    do q = 1, n_act_orb
!      qorb = list_act(q)
!      do p = 1, n_act_orb
!       porb = list_act(p)
       accu2 = 0d0
       do r = 1, n_act_orb
        rorb = list_act(r)
        do s = 1, n_act_orb
         sorb = list_act(s)
           do t= 1, n_act_orb
           torb = list_act(t)
           
           accu2 = accu2 &
                   + get_two_e_integral(porb,torb,rorb,sorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(r,s,q,t,istate) &
                   - get_two_e_integral(rorb,sorb,qorb,torb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(p,t,r,s,istate)


!  do q=1,mo_num ! list_act .....
!    do p=1,mo_num ! 
!      accu2 = 0d0
!      do r=1,mo_num
!        do s=1,mo_num
!          do t=1,mo_num
!            accu2 = accu2  &
!                      + get_two_e_integral(p,t,r,s,1) * y_act_2_rdm_spin_trace_mo(r,s,q,t,1) &
!                      - get_two_e_integral(r,s,q,t,1) * y_act_2_rdm_spin_trace_mo(p,t,r,s,1) 
                     ! - get_two_e_integral(q,t,r,s,1) *2d0 *act_2_rdm_spin_trace_mo(r,s,p,t,1) &
                      !+ get_two_e_integral(r,s,p,t,1) *2d0 *act_2_rdm_spin_trace_mo(q,t,r,s,1)
!          enddo
!        enddo
!      enddo
!      grad(p,q) = grad(p,q) + accu2
!    enddo
!  enddo

           enddo
          enddo
        enddo
        grad(p,q) = grad(p,q) + accu2 
      enddo
    enddo
  !enddo
 
! do q=1,mo_num
!   do p=1,q
!    grad(p,q) = -1d0 * (grad(p,q) - grad(q,p))
!    if (ABS(grad(p,q)) < 10d0**(-15)) then
!       grad(p,q) = 0d0
!     endif
!   enddo
! enddo
 
  ! Conversion mo_num*mo_num matrix to mo_num(mo_num-1)/2 vector
  !print*,'Gradient matrix -> vector'
  i=0
  do q = 1, mo_num
    do p = 1, q-1
      i=i+1
      v_grad(i) = -(grad(p,q) - grad(q,p))
    enddo
  enddo
  
  ! Debug
  !print*,'grad'
  !print*, grad(:,:)
  !print*, 'v_grad'
  !print*, v_grad(:)

  !==============
  ! Deallocation
  !==============

  deallocate(grad)

end subroutine
