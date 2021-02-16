subroutine hess(H,n)
  implicit none

  integer, intent(in) :: n 
  double precision, allocatable :: hessian(:,:,:,:)
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

  !do istate = 1, n_act_orb
  istate = 1  
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
                    hessian(p,q,r,s) = hessian(p,q,r,s) &
                                        + mo_one_e_integrals(porb,uorb) &
                                        * (one_e_dm_mo_alpha(u,s,istate) + one_e_dm_mo_beta(u,s,istate)) &
                                        + mo_one_e_integrals(uorb,sorb) &
                                        * (one_e_dm_mo_alpha(p,u,istate) + one_e_dm_mo_beta(p,u,istate))
                enddo
            endif
  
            if (p==s) then
                do u = 1, n_act_orb
                uorb = list_act(u)
                    hessian(p,q,r,s) = hessian(p,q,r,s) &
                                        + mo_one_e_integrals(uorb,rorb) &
                                        * (one_e_dm_mo_alpha(u,q,istate) + one_e_dm_mo_beta(u,q,istate)) &
                                        + mo_one_e_integrals(uorb,qorb) &
                                        * (one_e_dm_mo_alpha(r,u,istate) + one_e_dm_mo_beta(r,u,istate))
                enddo
            endif 
             
            hessian(p,q,r,s) = hessian(p,q,r,s) &
                                - mo_one_e_integrals(porb,sorb) &
                                * (one_e_dm_mo_alpha(r,q,istate) + one_e_dm_mo_beta(r,q,istate)) &
                                - mo_one_e_integrals(rorb,qorb) &
                                * (one_e_dm_mo_alpha(p,s,istate) + one_e_dm_mo_beta(p,s,istate))
  
  
            if (q==r) then
                do t = 1, n_act_orb
                torb = list_act(t)
                    do u = 1, n_act_orb
                    uorb = list_act(u)
                        do v = 1, n_act_orb
                        vorb = list_act(v)
                            hessian(p,q,r,s) = hessian(p,q,r,s) &
                             + get_two_e_integral(porb,torb,uorb,vorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(u,v,s,t,istate) &
                             + get_two_e_integral(uorb,vorb,sorb,torb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(p,t,u,v,istate)
                        enddo
                    enddo
                enddo
            endif
  
            if (p==s) then
                do t = 1, n_act_orb
                torb = list_act(t)
                    do u = 1, n_act_orb
                    uorb = list_act(u)
                        do v = 1, n_act_orb
                        vorb = list_act(vorb)
                            hessian(p,q,r,s) = hessian(p,q,r,s) &
                             + get_two_e_integral(uorb,vorb,qorb,torb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(r,t,u,v,istate) &
                             + get_two_e_integral(rorb,torb,uorb,vorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(u,v,q,t,istate)
                        enddo
                    enddo
                enddo
            endif
  
            do u = 1, n_act_orb
            uorb = list_act(u)
                do v = 1, n_act_orb
                vorb = list_act(v)
                    hessian(p,q,r,s) = hessian(p,q,r,s) &
                                        + get_two_e_integral(porb,rorb,uorb,vorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(u,v,q,s,istate) &
                                        + get_two_e_integral(uorb,vorb,qorb,sorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(p,r,u,v,istate)
                enddo
            enddo
  
            do t = 1, n_act_orb
            torb = list_act(t)
                do u = 1, n_act_orb
                uorb = list_act(u)
                    hessian(p,q,r,s) = hessian(p,q,r,s) &
                                        - get_two_e_integral(porb,uorb,sorb,torb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(r,t,q,u,istate) &
                                        - get_two_e_integral(porb,uorb,torb,sorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(t,r,q,u,istate) &
                                        - get_two_e_integral(rorb,torb,qorb,uorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(p,u,s,t,istate) &
                                        - get_two_e_integral(torb,rorb,qorb,uorb,mo_integrals_map) * y_act_2_rdm_spin_trace_mo(p,u,t,s,istate)
                enddo
            enddo
          enddo
        enddo
      enddo
    enddo 
  !enddo

  ! Debug
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
  do r = 1, mo_num
    do s = r+1, mo_num
      rs=rs+1
      pq=0
      do p = 1, mo_num 
        do q = p+1, mo_num  
            pq=pq+1
            print*,p,q,r,s
            H(pq,rs) = hessian(p,q,r,s) - hessian(q,p,r,s) - hessian(p,q,s,r) + hessian(q,p,s,r) ! pqrs-qprs-pqsr+qpsr
        enddo
      enddo
    enddo
  enddo

  ! Debug
  !print*,'H_pq,rs'
  !print*,H(:,:)
  !print*,'end'

  !==============
  ! Deallocation
  !==============

  deallocate(hessian)
  
end subroutine
