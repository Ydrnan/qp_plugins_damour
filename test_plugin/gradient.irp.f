subroutine gradient(grad)
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  
  ! Check if read_wf = true, else :
  ! qp set determinant read_wf true 
 
  END_DOC
  !double precision, allocatable :: grad(:,:)
  !double precision, allocatable :: v_grad(:)
  double precision, intent(out) :: grad(mo_num,mo_num)

  double precision :: get_two_e_integral
  double precision :: accu1, accu2
  integer :: i,n,p,q,r,s,t

  print*, mo_num
  
  !n=mo_num*(mo_num-1)/2

  !allocate(grad(mo_num,mo_num))
  !allocate(v_grad(n))

  grad = 0d0
  !v_grad = 0d0 
 
  do q=1,mo_num
    do p=1,mo_num
      accu1 = 0d0
      do r=1,mo_num
        accu1 = accu1 + mo_one_e_integrals(p,r) &
                       * (one_e_dm_mo_alpha(r,q,1) + one_e_dm_mo_beta(r,q,1)) &
                      - mo_one_e_integrals(r,q) &
                       * (one_e_dm_mo_alpha(p,r,1) + one_e_dm_mo_beta(p,r,1)) &
                      - mo_one_e_integrals(q,r) &
                       * (one_e_dm_mo_alpha(r,p,1) + one_e_dm_mo_beta(r,p,1)) &
                      + mo_one_e_integrals(r,p) &
                       * (one_e_dm_mo_alpha(q,r,1) + one_e_dm_mo_beta(q,r,1))
      enddo
      grad(p,q) = grad(p,q) +  accu1
    enddo
  enddo  


  do q=1,mo_num
    do p=1,mo_num
      accu2 = 0d0
      do r=1,mo_num
        do s=1,mo_num
          do t=1,mo_num
            accu2 = accu2  &
                      + 2d0 * get_two_e_integral(p,t,r,s,1) *act_2_rdm_spin_trace_mo(r,s,q,t,1) &
                      - 2d0 * get_two_e_integral(r,s,q,t,1) *act_2_rdm_spin_trace_mo(p,t,r,s,1) &
                      - 2d0 * get_two_e_integral(q,t,r,s,1) *act_2_rdm_spin_trace_mo(r,s,p,t,1) &
                      + 2d0 * get_two_e_integral(r,s,p,t,1) *act_2_rdm_spin_trace_mo(q,t,r,s,1)
          enddo
        enddo
      enddo
      grad(p,q) = grad(p,q) + accu2
    enddo
  enddo

 ! i = 0
 ! do q=2,mo_num
 !   do p=1,mo_num
 !     if (p<q) then
 !       v_grad(i) = grad(p,q)
 !       print*,p,q,grad(p,q)
 !       i = i+1 
 !     endif
 !   enddo
 !enddo
 
 do q=1,mo_num
   do p=1,mo_num
    grad(p,q) = 0.1d0 * grad(p,q) 
    if (ABS(grad(p,q)) < 10d0**(-15)) then
       grad(p,q) = 0d0
     endif
   enddo
 enddo
 
  do p=1,mo_num
      print*, grad(p,:)
  enddo

 !mo_coef(1,1) = 1d0
 ! Save the modification of the MOs
 !call save_mos
 !print*,mo_coef(1,1)

end
