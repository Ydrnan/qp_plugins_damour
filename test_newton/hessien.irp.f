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

  allocate(hessian(mo_num,mo_num,mo_num,mo_num))
  
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
    
  rs=0
  do r = 1, mo_num - 1
    do s = r+1, mo_num
      rs=rs+1
      pq=0
      do p = 1, mo_num - 1 
        do q = p+1, mo_num  
            pq=pq+1
            !print*,p,q,r,s
            H(pq,rs) = (hessian(p,q,r,s) - hessian(q,p,r,s) - hessian(p,q,s,r) + hessian(q,p,s,r)) ! pqrs-qprs-pqsr+qpsr
            !print*,p,q,r,s,H(pq,rs)
        enddo
      enddo
    enddo
  enddo

  allocate(h_tmpr(mo_num,mo_num,mo_num,mo_num))
  do r = 1, mo_num
    do s = 1, mo_num
      do q = 1, mo_num
        do p = 1, mo_num

          h_tmpr(p,q,r,s) = (hessian(p,q,r,s) - hessian(q,p,r,s) - hessian(p,q,s,r) + hessian(q,p,s,r))
          !if (ABS(toto(p,q,r,s)) > 1.d-10) then 
          !  toto(p,q,r,s) = toto(p,q,r,s)
          !else 
          !  toto(p,q,r,s) = 0d0
          !endif
        enddo
      enddo
    enddo
  enddo
  
  !print*,'Hessian'
  !open(unit=10,file='Hessien_test2_ifort.dat')
  !do p=1,mo_num
  !  do q=1,mo_num
  !     write(10,*) toto(p,q,:,:)
  !    print*, toto(p,q,:,:)
  !  enddo
  !enddo
  !close(10)
  !open(unit=10,file='one_body_dm_gfortran.dat')
  !do u=1,mo_num
  !do s=1,mo_num
  !write(10,*) (one_e_dm_mo_alpha(u,s,1) + one_e_dm_mo_beta(u,s,1))
  !enddo
  !enddo
  !close(10)

  !open(unit=11,file='two_body_dm_gfortran.dat')
  !do u=1,mo_num
  !do s=1,mo_num
  !write(11,*) two_e_dm_mo(u,s,:,:,1)
  !enddo
  !enddo
  !close(11)
  
  ! Debug
  !print*,'H_pq,rs'
  !print*,H(:,:)
  !print*,'end'

  !==============
  ! Deallocation
  !==============

  deallocate(hessian,h_tmpr)
  
end subroutine
