subroutine ost2aiai2_opt(nO,nV,det,t1_A,t2_A,M2_A)

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

  !! Deltas:((na, aa), (nb, jb))
    do bb = i_bb, f_bb
      if (bb == mb) cycle
      ba = bb - cc_nVa
  do ia = i_ia, f_ia
    if (ia == ma) cycle
    ib = ia + cc_nOa
      M2_A(ia,nb,na,bb) = M2_A(ia,nb,na,bb) &
      -1.0d0 * t2_A(ma, ib, ba, mb) &
      -1.0d0 * t1_A(ma, ba) * t1_A(ib, mb)
    enddo
  enddo

  !! Deltas:((mb, bb), (nb, jb))
    do aa = i_aa, f_aa
      if (aa == na) cycle
      ab = aa + cc_nVa
  do ia = i_ia, f_ia
    if (ia == ma) cycle
    ib = ia + cc_nOa
      M2_A(ia,nb,aa,mb) = M2_A(ia,nb,aa,mb) &
      -1.0d0 * t2_A(nb, ib, ab, mb) &
      -1.0d0 * t1_A(nb, ab) * t1_A(ib, mb)
    enddo
  enddo

  !! Deltas:((na, aa), (ma, ia))
    do bb = i_bb, f_bb
      if (bb == mb) cycle
      ba = bb - cc_nVa
  do jb = i_jb, f_jb
    if (jb == nb) cycle
    ja = jb - cc_nOa
      M2_A(ma,jb,na,bb) = M2_A(ma,jb,na,bb) &
      -1.0d0 * t2_A(ma, ja, ba, na) &
      -1.0d0 * t1_A(ma, ba) * t1_A(ja, na)
    enddo
  enddo

  !! Deltas:((mb, bb), (ma, ia))
    do aa = i_aa, f_aa
      if (aa == na) cycle
      ab = aa + cc_nVa
  do jb = i_jb, f_jb
    if (jb == nb) cycle
    ja = jb - cc_nOa
      M2_A(ma,jb,aa,mb) = M2_A(ma,jb,aa,mb) &
      -1.0d0 * t2_A(ja, nb, na, ab) &
      -1.0d0 * t1_A(nb, ab) * t1_A(ja, na)
    enddo
  enddo

  ! ### Spin case: i_a, j_b, a_b, b_a ###
  
  !! Deltas:((mb, ab), (nb, jb))
    do ba = i_ba, f_ba
      if (ba == na) cycle
      bb = ba + cc_nVa
  do ia = i_ia, f_ia
    if (ia == ma) cycle
    ib = ia + cc_nOa
      M2_A(ia,nb,mb,ba) = M2_A(ia,nb,mb,ba) &
      +1.0d0 * t2_A(nb, ib, bb, mb) &
      +1.0d0 * t1_A(nb, bb) * t1_A(ib, mb)
    enddo
  enddo

  !! Deltas:((na, ba), (nb, jb))
    do ab = i_ab, f_ab
      if (ab == mb) cycle
      aa = ab - cc_nVa
  do ia = i_ia, f_ia
    if (ia == ma) cycle
    ib = ia + cc_nOa
      M2_A(ia,nb,ab,na) = M2_A(ia,nb,ab,na) &
      +1.0d0 * t2_A(ma, ib, aa, mb) &
      +1.0d0 * t1_A(ma, aa) * t1_A(ib, mb)
    enddo
  enddo

  !! Deltas:((mb, ab), (ma, ia))
    do ba = i_ba, f_ba
      if (ba == na) cycle
      bb = ba + cc_nVa
  do jb = i_jb, f_jb
    if (jb == nb) cycle
    ja = jb - cc_nOa
      M2_A(ma,jb,mb,ba) = M2_A(ma,jb,mb,ba) &
      +1.0d0 * t2_A(ja, nb, na, bb) &
      +1.0d0 * t1_A(nb, bb) * t1_A(ja, na)
    enddo
  enddo

  !! Deltas:((na, ba), (ma, ia))
    do ab = i_ab, f_ab
      if (ab == mb) cycle
      aa = ab - cc_nVa
  do jb = i_jb, f_jb
    if (jb == nb) cycle
    ja = jb - cc_nOa
      M2_A(ma,jb,ab,na) = M2_A(ma,jb,ab,na) &
      +1.0d0 * t2_A(ma, ja, aa, na) &
      +1.0d0 * t1_A(ma, aa) * t1_A(ja, na)
    enddo
  enddo

  ! ### Spin case: i_b, j_a, a_a, b_b ###

  !! Deltas:((na, aa), (ma, ja))
    do bb = i_bb, f_bb
      if (bb == mb) cycle
      ba = bb - cc_nVa
  do ib = i_ib, f_ib
    if (ib == nb) cycle
    ia = ib - cc_nOa
      M2_A(ib,ma,na,bb) = M2_A(ib,ma,na,bb) &
      +1.0d0 * t2_A(ma, ia, ba, na) &
      +1.0d0 * t1_A(ma, ba) * t1_A(ia, na)
    enddo
  enddo

  !! Deltas:((mb, bb), (ma, ja))
    do aa = i_aa, f_aa
      if (aa == na) cycle
      ab = aa + cc_nVa
  do ib = i_ib, f_ib
    if (ib == nb) cycle
    ia = ib - cc_nOa
      M2_A(ib,ma,aa,mb) = M2_A(ib,ma,aa,mb) &
      +1.0d0 * t2_A(ia, nb, na, ab) &
      +1.0d0 * t1_A(nb, ab) * t1_A(ia, na)
    enddo
  enddo

  !! Deltas:((na, aa), (nb, ib))
    do bb = i_bb, f_bb
      if (bb == mb) cycle
      ba = bb - cc_nVa
  do ja = i_ja, f_ja
    if (ja == ma) cycle
    jb = ja + cc_nOa
      M2_A(nb,ja,na,bb) = M2_A(nb,ja,na,bb) &
      +1.0d0 * t2_A(ma, jb, ba, mb) &
      +1.0d0 * t1_A(ma, ba) * t1_A(jb, mb)
    enddo
  enddo

  !! Deltas:((mb, bb), (nb, ib))
    do aa = i_aa, f_aa
      if (aa == na) cycle
      ab = aa + cc_nVa
  do ja = i_ja, f_ja
    if (ja == ma) cycle
    jb = ja + cc_nOa
      M2_A(nb,ja,aa,mb) = M2_A(nb,ja,aa,mb) &
      +1.0d0 * t2_A(nb, jb, ab, mb) &
      +1.0d0 * t1_A(nb, ab) * t1_A(jb, mb)
    enddo
  enddo

  ! ### Spin case: i_b, j_a, a_b, b_a ###

  !! Deltas:((mb, ab), (ma, ja))
    do ba = i_ba, f_ba
      if (ba == na) cycle
      bb = ba + cc_nVa
  do ib = i_ib, f_ib
    if (ib == nb) cycle
    ia = ib - cc_nOa
      M2_A(ib,ma,mb,ba) = M2_A(ib,ma,mb,ba) &
      -1.0d0 * t2_A(ia, nb, na, bb) &
      -1.0d0 * t1_A(nb, bb) * t1_A(ia, na)
    enddo
  enddo

  !! Deltas:((na, ba), (ma, ja))
    do ab = i_ab, f_ab
      if (ab == mb) cycle
      aa = ab - cc_nVa
  do ib = i_ib, f_ib
    if (ib == nb) cycle
    ia = ib - cc_nOa
      M2_A(ib,ma,ab,na) = M2_A(ib,ma,ab,na) &
      -1.0d0 * t2_A(ma, ia, aa, na) &
      -1.0d0 * t1_A(ma, aa) * t1_A(ia, na)
    enddo
  enddo

  !! Deltas:((mb, ab), (nb, ib))
    do ba = i_ba, f_ba
      if (ba == na) cycle
      bb = ba + cc_nVa
  do ja = i_ja, f_ja
    if (ja == ma) cycle
    jb = ja + cc_nOa
      M2_A(nb,ja,mb,ba) = M2_A(nb,ja,mb,ba) &
      -1.0d0 * t2_A(nb, jb, bb, mb) &
      -1.0d0 * t1_A(nb, bb) * t1_A(jb, mb)
    enddo
  enddo

  !! Deltas:((na, ba), (nb, ib))
    do ab = i_ab, f_ab
      if (ab == mb) cycle
      aa = ab - cc_nVa
  do ja = i_ja, f_ja
    if (ja == ma) cycle
    jb = ja + cc_nOa
      M2_A(nb,ja,ab,na) = M2_A(nb,ja,ab,na) &
      -1.0d0 * t2_A(ma, jb, aa, mb) &
      -1.0d0 * t1_A(ma, aa) * t1_A(jb, mb)
    enddo
  enddo

end
