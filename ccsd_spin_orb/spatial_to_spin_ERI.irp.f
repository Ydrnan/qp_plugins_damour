BEGIN_PROVIDER [ double precision, dbERI, (spin_mo_num,spin_mo_num,spin_mo_num,spin_mo_num) ]
 implicit none
 BEGIN_DOC
 ! Anti-symmetrized Electron repulsion integrals in spin-orbital basis
 END_DOC

! Local variables

  integer                       :: p,q,r,s
  double precision,external     :: Kronecker_delta

! Output variables

  double precision, external :: get_two_e_integral

  double precision :: pqrs, pqsr

  PROVIDE mo_two_e_integrals_in_map

    do s=1,spin_mo_num
      do r=1,spin_mo_num
        do q=1,spin_mo_num
          do p=1,spin_mo_num
            pqrs = 0.d0
            pqsr = 0.d0
            if ( ( iand(p,1) == iand(r,1) )  .and. &
                 ( iand(q,1) == iand(s,1) ) ) then
!            pqrs = Kronecker_delta(mod(p,2),mod(r,2)) &
!                 * Kronecker_delta(mod(q,2),mod(s,2)) &
!                 * get_two_e_integral(                &
              pqrs = get_two_e_integral(                &
                      (p+1)/2,                         &
                      (q+1)/2,                         &
                      (r+1)/2,                         &
                      (s+1)/2,                         &
                      mo_two_e_integrals_in_map)
            endif

            if ( ( iand(p,1) == iand(s,1) )  .and. &
                 ( iand(q,1) == iand(r,1) ) ) then

!              pqsr = Kronecker_delta(mod(p,2),mod(s,2)) &
!                   * Kronecker_delta(mod(q,2),mod(r,2)) &
!                   * get_two_e_integral(                &
                pqsr = get_two_e_integral(                &
                      (p+1)/2,                         &
                      (q+1)/2,                         &
                      (s+1)/2,                         &
                      (r+1)/2,                         &
                      mo_two_e_integrals_in_map)
            endif

            dbERI(p,q,r,s) = pqrs - pqsr
          enddo
        enddo
      enddo
    enddo


END_PROVIDER


