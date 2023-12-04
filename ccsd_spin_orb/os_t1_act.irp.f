subroutine os_t1_act(nO,nV,det,t1_A,t2_A,M1_A)

  implicit none

  integer, intent(in)           :: nO,nV
  integer(bit_kind), intent(in) :: det(N_int,2)
  double precision, intent(in)  :: t1_A(nO,nV), t2_A(nO,nO,nV,nV)
  
  double precision, intent(out) :: M1_A(nO,nV)

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

  ! ### Spin case: i_a, a_a ###

  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    do ia = i_ia, f_ia
      if (ia == ma) cycle 
      ib = ia + cc_nOa
      M1_A(ia,aa) = M1_A(ia,aa) & 
      -1.0d0 * t1_A(nb, ab) * t2_A(ma, ib, na, mb) & 
      -1.0d0 * t1_A(ma, na) * t2_A(nb, ib, ab, mb) & 
      -1.0d0 * t1_A(ib, mb) * t2_A(ma, nb, na, ab) & 
      +1.0d0 * t1_A(nb, mb) * t2_A(ma, ib, na, ab) & 
      -1.0d0 * t1_A(nb, ab) * t1_A(ma, na) * t1_A(ib, mb)
    enddo
  enddo

  !! Deltas:((na, aa))
  do ia = i_ia, f_ia
    if (ia == ma) cycle 
    ib = ia + cc_nOa
    M1_A(ia,na) = M1_A(ia,na) & 
    -1.0d0 * t2_A(ma, ib, na, mb) & 
    -1.0d0 * t1_A(ma, na) * t1_A(ib, mb)
  enddo

  !! Deltas:((ma, ia))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    M1_A(ma,aa) = M1_A(ma,aa) & 
    +1.0d0 * t2_A(ma, nb, na, ab) & 
    +1.0d0 * t1_A(nb, ab) * t1_A(ma, na)
  enddo

  !! Deltas:((ma, ia), (na, aa))
  M1_A(ma,na) = M1_A(ma,na) & 
  +1.0d0 * t1_A(ma, na)

  ! ### Spin case: i_b, a_b ###

  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    do ib = i_ib, f_ib
      if (ib == nb) cycle 
      ia = ib - cc_nOa
      M1_A(ib,ab) = M1_A(ib,ab) & 
      -1.0d0 * t1_A(ma, aa) * t2_A(ia, nb, na, mb) & 
      -1.0d0 * t1_A(ia, na) * t2_A(ma, nb, aa, mb) & 
      +1.0d0 * t1_A(ma, na) * t2_A(ia, nb, aa, mb) & 
      -1.0d0 * t1_A(nb, mb) * t2_A(ma, ia, aa, na) & 
      -1.0d0 * t1_A(ma, aa) * t1_A(ia, na) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((mb, ab))
  do ib = i_ib, f_ib
    if (ib == nb) cycle 
    ia = ib - cc_nOa
    M1_A(ib,mb) = M1_A(ib,mb) & 
    -1.0d0 * t2_A(ia, nb, na, mb) & 
    -1.0d0 * t1_A(ia, na) * t1_A(nb, mb)
  enddo

  !! Deltas:((nb, ib))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    M1_A(nb,ab) = M1_A(nb,ab) & 
    +1.0d0 * t2_A(ma, nb, aa, mb) & 
    +1.0d0 * t1_A(ma, aa) * t1_A(nb, mb)
  enddo

  !! Deltas:((nb, ib), (mb, ab))
  M1_A(nb,mb) = M1_A(nb,mb) & 
  +1.0d0 * t1_A(nb, mb)

end
