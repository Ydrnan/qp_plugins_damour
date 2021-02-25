subroutine hess(n,H)
  implicit none

  !==================================================================
  ! Compute the hessian of energy with respects to orbital rotations
  !==================================================================

  integer, intent(in) :: n ! mo_num*(mo_num-1)/2
  double precision, allocatable :: hessian(:,:,:,:), h_tmpr(:,:,:,:)
  double precision, intent(out) :: H(n,n)
  integer :: p,q
  integer :: r,s,t,u,v
  integer :: pq,rs
  integer :: istate
  double precision :: get_two_e_integral

  ! Provided :
  ! mo_one_e_integrals : mono e- integrals
  ! get_two_e_integral : two e- integrals
  ! one_e_dm_mo_alpha, one_e_dm_mo_beta : one body density matrix
  ! two_e_dm_mo : two body density matrix

  !============
  ! Allocation
  !============

  allocate(hessian(mo_num,mo_num,mo_num,mo_num),h_tmpr(mo_num,mo_num,mo_num,mo_num))

  !=============
  ! Calculation
  !=============

  hessian = 0d0

  !do istate = 1, N_states
  istate = 1

! First line, first term
    do p = 1, mo_num
      do q = 1, mo_num
        do r = 1, mo_num
          do s = 1, mo_num

              if (q==r) then
                do u = 1, mo_num

                  hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * (  &
                    mo_one_e_integrals(u,p) * (one_e_dm_mo_alpha(u,s,istate) + one_e_dm_mo_beta(u,s,istate)) &
                  + mo_one_e_integrals(s,u) * (one_e_dm_mo_alpha(p,u,istate) + one_e_dm_mo_beta(p,u,istate)))

                enddo
              endif

          enddo
        enddo
      enddo
    enddo

 ! First line, second term
     do p = 1, mo_num
       do q = 1, mo_num
         do r = 1, mo_num
           do s = 1, mo_num

             if (p==s) then
               do u = 1, mo_num

                     hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * ( &
                       mo_one_e_integrals(u,r) * (one_e_dm_mo_alpha(u,q,istate) + one_e_dm_mo_beta(u,q,istate)) &
                     + mo_one_e_integrals(q,u) * (one_e_dm_mo_alpha(r,u,istate) + one_e_dm_mo_beta(r,u,istate)))
               enddo
             endif

           enddo
         enddo
       enddo
     enddo

 !! First line, third term
     do p = 1, mo_num
       do q = 1, mo_num
         do r = 1, mo_num
           do s = 1, mo_num

             hessian(p,q,r,s) = hessian(p,q,r,s) &
             - mo_one_e_integrals(s,p) * (one_e_dm_mo_alpha(r,q,istate) + one_e_dm_mo_beta(r,q,istate)) &
             - mo_one_e_integrals(q,r) * (one_e_dm_mo_alpha(p,s,istate) + one_e_dm_mo_beta(p,s,istate))

           enddo
         enddo
       enddo
     enddo

 !Second line, first term
     do p = 1, mo_num
       do q = 1, mo_num
         do r = 1, mo_num
           do s = 1, mo_num

              if (q==r) then
                do t = 1, mo_num
                  do u = 1, mo_num
                    do v = 1, mo_num

                      hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * (  &
                        get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,s,t,istate) &
                      + get_two_e_integral(s,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v,istate))

                    enddo
                  enddo
                enddo
              endif

           enddo
         enddo
       enddo
     enddo

! Second line, second term
     do p = 1, mo_num
       do q = 1, mo_num
         do r = 1, mo_num
           do s = 1, mo_num

             if (p==s) then
               do t = 1, mo_num
                 do u = 1, mo_num
                   do v = 1, mo_num

                     hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * ( &
                       get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(r,t,u,v,istate) &
                     + get_two_e_integral(u,v,r,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t,istate))

                   enddo
                 enddo
               enddo
             endif

           enddo
         enddo
       enddo
     enddo

! Third line, first term
    do p = 1, mo_num
      do q = 1, mo_num
        do r = 1, mo_num
          do s = 1, mo_num

            do u = 1, mo_num
              do v = 1, mo_num

                hessian(p,q,r,s) = hessian(p,q,r,s) &
                 + get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s,istate) &
                 + get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v,istate)

              enddo
            enddo

          enddo
        enddo
      enddo
    enddo

! Third line, second term
    do p = 1, mo_num
      do q = 1, mo_num
        do r = 1, mo_num
          do s = 1, mo_num

            do t = 1, mo_num
              do u = 1, mo_num

                hessian(p,q,r,s) = hessian(p,q,r,s) &
                 - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u,istate) &
                 - get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u,istate) &
                 - get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t,istate) &
                 - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s,istate)

              enddo
            enddo

          enddo
        enddo
      enddo
    enddo

  !enddo

  !===========
  ! 2D matrix
  !===========

  ! Convert the hessian mo_num * mo_num * mo_num * mo_num matrix in a
  ! 2D n * n matrix (n = mo_num*(mo_num-1)/2)
  ! H(pq,rs) : p<q and r<s
  ! Hessian(p,q,r,s) = P_pq P_rs [ ...]
  ! => Hessian(p,q,r,s) = (p,q,r,s) - (q,p,r,s) - (p,q,s,r) + (q,p,s,r)

!  rs=0
!  do r = 1, mo_num
!    do s = 1, r-1
!      rs=rs+1
!      pq=0
!      do p = 1, mo_num
!        do q = 1, p-1
!            pq=pq+1
!
!            H(pq,rs) = (hessian(p,q,r,s) - hessian(q,p,r,s) - hessian(p,q,s,r) + hessian(q,p,s,r)) ! pqrs-qprs-pqsr+qpsr
!            H(pq,rs)=0.5d0*H(pq,rs)
!            ! TODO : check why 2.0 is here
!
!        enddo
!      enddo
!    enddo
!  enddo

  ! Verification des éléments du hessien
  do r = 1, mo_num
    do s = 1, mo_num
      do q = 1, mo_num
        do p = 1, mo_num

          !h_tmpr(p,q,r,s) =0.5d0* (hessian(p,q,r,s) - hessian(q,p,r,s) - hessian(p,q,s,r) + hessian(q,p,s,r))
          h_tmpr(p,q,r,s) = (hessian(p,q,r,s) - hessian(q,p,r,s) - hessian(p,q,s,r) + hessian(q,p,s,r))

        enddo
      enddo
    enddo
  enddo

  double precision :: e_val(mo_num**2),H_v(mo_num**2,mo_num**2), H_u(mo_num,mo_num,mo_num,mo_num)
  call lapack_diag(e_val,H_v,h_tmpr,mo_num**2,mo_num**2)

  print*,'e_val'
  write(*,'(100(F10.5))') e_val(:) 

!  print*,'verif H'
!  do p = 1, mo_num
!    do q = 1, mo_num
!          write(*,'(100(F10.5))') h_tmpr(p,q,1:mo_num,1:mo_num)
!    enddo
!  enddo

  do pq = 1, n
    call in_mat_vec_index(pq,p,q)
    do rs = 1, n
      call in_mat_vec_index(rs,r,s)
      H(pq,rs) = h_tmpr(p,q,r,s)   
    enddo
  enddo

  !==============
  ! Deallocation
  !==============

  deallocate(hessian,h_tmpr)

end subroutine
