subroutine hess(n,H)
  implicit none

  integer, intent(in) :: n 
  double precision, allocatable :: hessian(:,:,:,:), toto(:,:,:,:)
  double precision, intent(out) :: H(n,n)
  integer :: p,q
  integer :: r,s,t,u,v
  integer :: pq,rs
  integer :: porb,qorb,rorb,sorb,torb,uorb,vorb
  integer :: istate
  double precision :: get_two_e_integral

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
!    do p = 1, n_act_orb
!    porb = list_act(p) 
!      do q = 1, n_act_orb
!      qorb = list_act(q)
!        do r = 1, n_act_orb
!        rorb = list_act(r)  
!          do s = 1, n_act_orb
!          sorb = list_act(s)

! First line, first term 
    do p = 1, n_act_orb
    porb = list_act(p)
      do q = 1, n_act_orb
      qorb = list_act(q)
        do r = 1, n_act_orb
        rorb = list_act(r)
          do s = 1, n_act_orb
          sorb = list_act(s) 

              if (q==r) then
                do u = 1, n_act_orb
                uorb = list_act(u)

                  hessian(p,q,q,s) = hessian(p,q,q,s) + 0.5d0 * (  &
                    mo_one_e_integrals(uorb,porb) * (one_e_dm_mo_alpha(u,s,istate) + one_e_dm_mo_beta(u,s,istate)) &
                  + mo_one_e_integrals(sorb,uorb) * (one_e_dm_mo_alpha(p,u,istate) + one_e_dm_mo_beta(p,u,istate)))

                enddo
              endif

          enddo
        enddo
      enddo
    enddo

 ! First line, second term
     do p = 1, n_act_orb
     porb = list_act(p)
       do q = 1, n_act_orb
       qorb = list_act(q)
         do r = 1, n_act_orb
         rorb = list_act(r)
           do s = 1, n_act_orb
           sorb = list_act(s)
 
             if (p==s) then
               do u = 1, n_act_orb
                 uorb = list_act(u)
                     hessian(p,q,r,p) = hessian(p,q,r,p) + 0.5d0 * ( &
                       mo_one_e_integrals(uorb,rorb) * (one_e_dm_mo_alpha(u,q,istate) + one_e_dm_mo_beta(u,q,istate)) &
                     + mo_one_e_integrals(qorb,uorb) * (one_e_dm_mo_alpha(r,u,istate) + one_e_dm_mo_beta(r,u,istate)))
               enddo
             endif
 
           enddo
         enddo
       enddo
     enddo
 
 !! First line, third term
     do p = 1, n_act_orb
     porb = list_act(p)
       do q = 1, n_act_orb
       qorb = list_act(q)
         do r = 1, n_act_orb
         rorb = list_act(r)
           do s = 1, n_act_orb
           sorb = list_act(s)
              
             hessian(p,q,r,s) = hessian(p,q,r,s) &
             - mo_one_e_integrals(sorb,porb) * (one_e_dm_mo_alpha(r,q,istate) + one_e_dm_mo_beta(r,q,istate)) &
             - mo_one_e_integrals(qorb,rorb) * (one_e_dm_mo_alpha(p,s,istate) + one_e_dm_mo_beta(p,s,istate))
 
           enddo
         enddo
       enddo
     enddo
 
 !Second line, first term
     do p = 1, n_act_orb
     porb = list_act(p)
       do q = 1, n_act_orb
       qorb = list_act(q)
         do r = 1, n_act_orb
         rorb = list_act(r)
           do s = 1, n_act_orb
           sorb = list_act(s)
             
              if (q==r) then
                do t = 1, n_act_orb
                torb = list_act(t)
                  do u = 1, n_act_orb
                  uorb = list_act(u)
                    do v = 1, n_act_orb
                    vorb =list_act(v)
 
                      hessian(p,q,q,s) = hessian(p,q,q,s) + 0.5d0 * (  &
                        get_two_e_integral(uorb,vorb,porb,torb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(u,v,s,t,istate) &
                      + get_two_e_integral(sorb,torb,uorb,vorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(p,t,u,v,istate))
 
                    enddo
                  enddo
                enddo
              endif
 
           enddo
         enddo
       enddo
     enddo

! Second line, second term
     do p = 1, n_act_orb
     porb = list_act(p)
       do q = 1, n_act_orb
       qorb = list_act(q)
         do r = 1, n_act_orb
         rorb = list_act(r)
           do s = 1, n_act_orb
           sorb = list_act(s)
 
             if (p==s) then
               do t = 1, n_act_orb
               torb = list_act(t)
                 do u = 1, n_act_orb
                 uorb = list_act(u)
                   do v = 1, n_act_orb
                   vorb = list_act(v)
 
                     hessian(p,q,r,p) = hessian(p,q,r,p) + 0.5d0 * ( &
                       get_two_e_integral(qorb,torb,uorb,vorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(r,t,u,v,istate) &
                     + get_two_e_integral(uorb,vorb,rorb,torb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(u,v,q,t,istate))
 
                   enddo
                 enddo
               enddo
             endif
 
           enddo
         enddo
       enddo
     enddo

! Third line, first term
    do p = 1, n_act_orb
    porb = list_act(p)
      do q = 1, n_act_orb
      qorb = list_act(q)
        do r = 1, n_act_orb
        rorb = list_act(r)
          do s = 1, n_act_orb
          sorb = list_act(s)

            do u = 1, n_act_orb
            uorb = list_act(u)
              do v = 1, n_act_orb
              vorb = list_act(v)

                hessian(p,q,r,s) = hessian(p,q,r,s) &
                 + get_two_e_integral(uorb,vorb,porb,rorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(u,v,q,s,istate) &
                 + get_two_e_integral(qorb,sorb,uorb,vorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(p,r,u,v,istate)

              enddo
            enddo

          enddo
        enddo
      enddo
    enddo

! Third line, second term
    do p = 1, n_act_orb
    porb = list_act(p)
      do q = 1, n_act_orb
      qorb = list_act(q)
        do r = 1, n_act_orb
        rorb = list_act(r)
          do s = 1, n_act_orb
          sorb = list_act(s)
  
            do t = 1, n_act_orb
            torb = list_act(t)
              do u = 1, n_act_orb
              uorb = list_act(u)

                hessian(p,q,r,s) = hessian(p,q,r,s) &
                 - get_two_e_integral(sorb,torb,porb,uorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(r,t,q,u,istate) &
                 - get_two_e_integral(torb,sorb,porb,uorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(t,r,q,u,istate) &
                 - get_two_e_integral(qorb,uorb,rorb,torb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(p,u,s,t,istate) &
                 - get_two_e_integral(qorb,uorb,torb,rorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(p,u,t,s,istate)

              enddo
            enddo

          enddo
        enddo
      enddo
    enddo
 
  !enddo

  ! Debug
  !print*,'two bd density'
  !print*,y_act_2_rdm_spin_trace_mo(:,:,:,:,1)
  ! Le hessien est focément faux vu que les permutations sont faites après 
  !print*,'Hessian'
  !print*,hessian(:,:,:,:)

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
            H(pq,rs) = hessian(p,q,r,s) - hessian(q,p,r,s) - hessian(p,q,s,r) + hessian(q,p,s,r) ! pqrs-qprs-pqsr+qpsr
            !print*,p,q,r,s,H(pq,rs)
        enddo
      enddo
    enddo
  enddo

  allocate(toto(mo_num,mo_num,mo_num,mo_num))
  do r = 1, mo_num
    do s = 1, mo_num
      do q = 1, mo_num
        do p = 1, mo_num

          toto(p,q,r,s) = 0.5d0*(hessian(p,q,r,s) - hessian(q,p,r,s) - hessian(p,q,s,r) + hessian(q,p,s,r))
          if (ABS(toto(p,q,r,s)) > 10d0**(-15)) then 
          toto(p,q,r,s) = toto(p,q,r,s)
          else 
          toto(p,q,r,s) = 0d0
          endif
        enddo
      enddo
    enddo
  enddo
  
  print*,'Hessian'
  open(unit=10,file='Hessien_test.dat')
  do p=1,mo_num
    do q=1,mo_num
      print*, toto(p,q,:,:)
    enddo
  enddo
  close(10)
  ! Debug
  !print*,'H_pq,rs'
  !print*,H(:,:)
  !print*,'end'

  !==============
  ! Deallocation
  !==============

  deallocate(hessian)
  
end subroutine
