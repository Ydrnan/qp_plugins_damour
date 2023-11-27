subroutine ost2aaii_opt(nO,nV,det,t1_A,t2_A,M2_A)

  implicit none

  integer, intent(in)           :: nO,nV
  integer(bit_kind), intent(in) :: det(N_int,2)
  double precision, intent(in)  :: t1_A(nO,nV), t2_A(nO,nO,nV,nV)

  double precision, intent(out) :: M2_A(nO,nO,nV,nV)

  integer                       :: ia,ib,ja,jb,na,nb,ma,mb,aa,ab,ba,bb
  integer                       :: i_ia, i_ja, i_aa, i_ba
  integer                       :: i_ib, i_jb, i_ab, i_bb
  integer                       :: f_ia, f_ja, f_aa, f_ba
  integer                       :: f_ib, f_jb, f_ab, f_bb

  ! List of open spin orbitals
  call extract_open_spin_orb(nO,nV,det,ma,mb,na,nb)

  i_ia = 1
  i_ja = 1
  i_ib = cc_nOa + 1
  i_jb = cc_nOa + 1
  i_aa = 1
  i_ba = 1
  i_ab = cc_nVa + 1
  i_bb = cc_nVa + 1

  f_ia = cc_nOa
  f_ja = cc_nOa
  f_ib = cc_nOab
  f_jb = cc_nOab
  f_aa = cc_nVa
  f_ba = cc_nVa
  f_ab = cc_nVab
  f_bb = cc_nVab

  ! ### Spin case: i_a, j_b, a_a, b_b ###
 
  !! Deltas:((ma, ia), (nb, jb))
    do bb = i_bb, f_bb
      if (bb == mb) cycle
      ba = bb - cc_nVa
  do aa = i_aa, f_aa
    if (aa == na) cycle
    ab = aa + cc_nVa
      M2_A(ma,nb,aa,bb) = M2_A(ma,nb,aa,bb) &
      +1.0d0 * t2_A(ma, nb, ba, ab) &
      +1.0d0 * t1_A(nb, ab) * t1_A(ma, ba)
    enddo
  enddo
  
  ! ### Spin case: i_a, j_b, a_b, b_a ###

  !! Deltas:((ma, ia), (nb, jb))
    do ba = i_ba, f_ba
      if (ba == na) cycle
      bb = ba + cc_nVa
  do ab = i_ab, f_ab
    if (ab == mb) cycle
    aa = ab - cc_nVa
      M2_A(ma,nb,ab,ba) = M2_A(ma,nb,ab,ba) &
      -1.0d0 * t2_A(ma, nb, aa, bb) &
      -1.0d0 * t1_A(ma, aa) * t1_A(nb, bb)
    enddo
  enddo

  ! ### Spin case: i_b, j_a, a_a, b_b ###

  !! Deltas:((ma, ja), (nb, ib))
    do bb = i_bb, f_bb
      if (bb == mb) cycle
      ba = bb - cc_nVa
  do aa = i_aa, f_aa
    if (aa == na) cycle
    ab = aa + cc_nVa
      M2_A(nb,ma,aa,bb) = M2_A(nb,ma,aa,bb) &
      -1.0d0 * t2_A(ma, nb, ba, ab) &
      -1.0d0 * t1_A(nb, ab) * t1_A(ma, ba)
    enddo
  enddo

  ! ### Spin case: i_b, j_a, a_b, b_a ###

  !! Deltas:((ma, ja), (nb, ib))
    do ba = i_ba, f_ba
      if (ba == na) cycle
      bb = ba + cc_nVa
  do ab = i_ab, f_ab
    if (ab == mb) cycle
    aa = ab - cc_nVa
      M2_A(nb,ma,ab,ba) = M2_A(nb,ma,ab,ba) &
      +1.0d0 * t2_A(ma, nb, aa, bb) &
      +1.0d0 * t1_A(ma, aa) * t1_A(nb, bb)
    enddo
  enddo
  


end
