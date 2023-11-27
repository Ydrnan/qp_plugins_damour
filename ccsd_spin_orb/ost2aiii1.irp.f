subroutine ost2aiii1(nO,nV,det,t1_A,t2_A,M2_A)

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

  ! ### Spin case: i_a, j_a, a_a, b_a ###

  !! Deltas:((ma, ja))
  do ia = i_ia, f_ia
    if (ia == ma) cycle
    ib = ia + cc_nOa
    do aa = i_aa, f_aa
      if (aa == na) cycle
      ab = aa + cc_nVa
      do ba = i_ba, f_ba
        if (ba == na) cycle
        bb = ba + cc_nVa
        M2_A(ia,ma,aa,ba) = M2_A(ia,ma,aa,ba) &
        -1.0d0 * t1_A(nb, ab) * t2_A(ma, ib, na, bb) &
        +1.0d0 * t1_A(nb, bb) * t2_A(ma, ib, na, ab)
      enddo
    enddo
  enddo

  !! Deltas:((ma, ia))
  do ja = i_ja, f_ja
    if (ja == ma) cycle
    jb = ja + cc_nOa
    do aa = i_aa, f_aa
      if (aa == na) cycle
      ab = aa + cc_nVa
      do ba = i_ba, f_ba
        if (ba == na) cycle
        bb = ba + cc_nVa
        M2_A(ma,ja,aa,ba) = M2_A(ma,ja,aa,ba) &
        +1.0d0 * t1_A(nb, ab) * t2_A(ma, jb, na, bb) &
        -1.0d0 * t1_A(nb, bb) * t2_A(ma, jb, na, ab)
      enddo
    enddo
  enddo

  ! ### Spin case: i_a, j_b, a_a, b_b ###
  !! Deltas:((nb, jb))
  do ia = i_ia, f_ia
    if (ia == ma) cycle
    ib = ia + cc_nOa
    do aa = i_aa, f_aa
      if (aa == na) cycle
      ab = aa + cc_nVa
      do bb = i_bb, f_bb
        if (bb == mb) cycle
        ba = bb - cc_nVa
        M2_A(ia,nb,aa,bb) = M2_A(ia,nb,aa,bb) &
        -1.0d0 * t1_A(nb, ab) * t2_A(ma, ib, ba, mb) &
        -1.0d0 * t1_A(ma, ba) * t2_A(nb, ib, ab, mb)
      enddo
    enddo
  enddo

  !! Deltas:((ma, ia))
  do jb = i_jb, f_jb
    if (jb == nb) cycle
    ja = jb - cc_nOa
    do aa = i_aa, f_aa
      if (aa == na) cycle
      ab = aa + cc_nVa
      do bb = i_bb, f_bb
        if (bb == mb) cycle
        ba = bb - cc_nVa
        M2_A(ma,jb,aa,bb) = M2_A(ma,jb,aa,bb) &
        -1.0d0 * t1_A(nb, ab) * t2_A(ma, ja, ba, na) &
        -1.0d0 * t1_A(ma, ba) * t2_A(ja, nb, na, ab)
      enddo
    enddo
  enddo

  ! ### Spin case: i_a, j_b, a_b, b_a ###
  !! Deltas:((nb, jb))
  do ia = i_ia, f_ia
    if (ia == ma) cycle
    ib = ia + cc_nOa
    do ab = i_ab, f_ab
      if (ab == mb) cycle
      aa = ab - cc_nVa
      do ba = i_ba, f_ba
        if (ba == na) cycle
        bb = ba + cc_nVa
        M2_A(ia,nb,ab,ba) = M2_A(ia,nb,ab,ba) &
        +1.0d0 * t1_A(ma, aa) * t2_A(nb, ib, bb, mb) &
        +1.0d0 * t1_A(nb, bb) * t2_A(ma, ib, aa, mb)
      enddo
    enddo
  enddo

  !! Deltas:((ma, ia))
  do jb = i_jb, f_jb
    if (jb == nb) cycle
    ja = jb - cc_nOa
    do ab = i_ab, f_ab
      if (ab == mb) cycle
      aa = ab - cc_nVa
      do ba = i_ba, f_ba
        if (ba == na) cycle
        bb = ba + cc_nVa
        M2_A(ma,jb,ab,ba) = M2_A(ma,jb,ab,ba) &
        +1.0d0 * t1_A(ma, aa) * t2_A(ja, nb, na, bb) &
        +1.0d0 * t1_A(nb, bb) * t2_A(ma, ja, aa, na)
      enddo
    enddo
  enddo

  ! ### Spin case: i_b, j_a, a_a, b_b ###
  !! Deltas:((ma, ja))
  do ib = i_ib, f_ib
    if (ib == nb) cycle
    ia = ib - cc_nOa
    do aa = i_aa, f_aa
      if (aa == na) cycle
      ab = aa + cc_nVa
      do bb = i_bb, f_bb
        if (bb == mb) cycle
        ba = bb - cc_nVa
        M2_A(ib,ma,aa,bb) = M2_A(ib,ma,aa,bb) &
        +1.0d0 * t1_A(nb, ab) * t2_A(ma, ia, ba, na) &
        +1.0d0 * t1_A(ma, ba) * t2_A(ia, nb, na, ab)
      enddo
    enddo
  enddo

  !! Deltas:((nb, ib))
  do ja = i_ja, f_ja
    if (ja == ma) cycle
    jb = ja + cc_nOa
    do aa = i_aa, f_aa
      if (aa == na) cycle
      ab = aa + cc_nVa
      do bb = i_bb, f_bb
        if (bb == mb) cycle
        ba = bb - cc_nVa
        M2_A(nb,ja,aa,bb) = M2_A(nb,ja,aa,bb) &
        +1.0d0 * t1_A(nb, ab) * t2_A(ma, jb, ba, mb) &
        +1.0d0 * t1_A(ma, ba) * t2_A(nb, jb, ab, mb)
      enddo
    enddo
  enddo

  ! ### Spin case: i_b, j_a, a_b, b_a ###
  !! Deltas:((ma, ja))
  do ib = i_ib, f_ib
    if (ib == nb) cycle
    ia = ib - cc_nOa
    do ab = i_ab, f_ab
      if (ab == mb) cycle
      aa = ab - cc_nVa
      do ba = i_ba, f_ba
        if (ba == na) cycle
        bb = ba + cc_nVa
        M2_A(ib,ma,ab,ba) = M2_A(ib,ma,ab,ba) &
        -1.0d0 * t1_A(ma, aa) * t2_A(ia, nb, na, bb) &
        -1.0d0 * t1_A(nb, bb) * t2_A(ma, ia, aa, na)
      enddo
    enddo
  enddo

  !! Deltas:((nb, ib))
  do ja = i_ja, f_ja
    if (ja == ma) cycle
    jb = ja + cc_nOa
    do ab = i_ab, f_ab
      if (ab == mb) cycle
      aa = ab - cc_nVa
      do ba = i_ba, f_ba
        if (ba == na) cycle
        bb = ba + cc_nVa
        M2_A(nb,ja,ab,ba) = M2_A(nb,ja,ab,ba) &
        -1.0d0 * t1_A(ma, aa) * t2_A(nb, jb, bb, mb) &
        -1.0d0 * t1_A(nb, bb) * t2_A(ma, jb, aa, mb)
      enddo
    enddo
  enddo

  ! ### Spin case: i_b, j_b, a_b, b_b ###
  !! Deltas:((nb, jb))
  do ib = i_ib, f_ib
    if (ib == nb) cycle
    ia = ib - cc_nOa
    do ab = i_ab, f_ab
      if (ab == mb) cycle
      aa = ab - cc_nVa
      do bb = i_bb, f_bb
        if (bb == mb) cycle
        ba = bb - cc_nVa
        M2_A(ib,nb,ab,bb) = M2_A(ib,nb,ab,bb) &
        -1.0d0 * t1_A(ma, aa) * t2_A(ia, nb, ba, mb) &
        +1.0d0 * t1_A(ma, ba) * t2_A(ia, nb, aa, mb)
      enddo
    enddo
  enddo

  !! Deltas:((nb, ib))
  do jb = i_jb, f_jb
    if (jb == nb) cycle
    ja = jb - cc_nOa
    do ab = i_ab, f_ab
      if (ab == mb) cycle
      aa = ab - cc_nVa
      do bb = i_bb, f_bb
        if (bb == mb) cycle
        ba = bb - cc_nVa
        M2_A(nb,jb,ab,bb) = M2_A(nb,jb,ab,bb) &
        +1.0d0 * t1_A(ma, aa) * t2_A(ja, nb, ba, mb) &
        -1.0d0 * t1_A(ma, ba) * t2_A(ja, nb, aa, mb)
      enddo
    enddo
  enddo

end
