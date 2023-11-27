subroutine ost2iiii1(nO,nV,det,t1_A,t2_A,M2_A)

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

  do ia = i_ia, f_ia
    if (ia == ma) cycle
    ib = ia + cc_nOa
    do ja = i_ja, f_ja
      if (ja == ma) cycle
      jb = ja + cc_nOa
      do aa = i_aa, f_aa
        if (aa == na) cycle
        ab = aa + cc_nVa
        do ba = i_ba, f_ba
          if (ba == na) cycle
          bb = ba + cc_nVa
          M2_A(ia,ja,aa,ba) = M2_A(ia,ja,aa,ba) &
          -1.0d0 * t2_A(nb, jb, ab, bb) * t2_A(ma, ib, na, mb) &
          +1.0d0 * t2_A(nb, ib, ab, bb) * t2_A(ma, jb, na, mb)
        enddo
      enddo
    enddo
  enddo

    ! ### Spin case: i_a, j_b, a_a, b_b ###

  do ia = i_ia, f_ia
    if (ia == ma) cycle
    ib = ia + cc_nOa
    do jb = i_jb, f_jb
      if (jb == nb) cycle
      ja = jb - cc_nOa
      do aa = i_aa, f_aa
        if (aa == na) cycle
        ab = aa + cc_nVa
        do bb = i_bb, f_bb
          if (bb == mb) cycle
          ba = bb - cc_nVa
          M2_A(ia,jb,aa,bb) = M2_A(ia,jb,aa,bb) &
          -1.0d0 * t2_A(ja, nb, ba, ab) * t2_A(ma, ib, na, mb) &
          -1.0d0 * t2_A(ma, ib, ba, ab) * t2_A(ja, nb, na, mb)
        enddo
      enddo
    enddo
  enddo


  ! ### Spin case: i_a, j_b, a_b, b_a ###

  do ia = i_ia, f_ia
    if (ia == ma) cycle
    ib = ia + cc_nOa
    do jb = i_jb, f_jb
      if (jb == nb) cycle
      ja = jb - cc_nOa
      do ab = i_ab, f_ab
        if (ab == mb) cycle
        aa = ab - cc_nVa
        do ba = i_ba, f_ba
          if (ba == na) cycle
          bb = ba + cc_nVa
          M2_A(ia,jb,ab,ba) = M2_A(ia,jb,ab,ba) &
          +1.0d0 * t2_A(ja, nb, aa, bb) * t2_A(ma, ib, na, mb) &
          +1.0d0 * t2_A(ma, ib, aa, bb) * t2_A(ja, nb, na, mb)
        enddo
      enddo
    enddo
  enddo

  ! ### Spin case: i_b, j_a, a_a, b_b ###

  do ib = i_ib, f_ib
    if (ib == nb) cycle
    ia = ib - cc_nOa
    do ja = i_ja, f_ja
      if (ja == ma) cycle
      jb = ja + cc_nOa
      do aa = i_aa, f_aa
        if (aa == na) cycle
        ab = aa + cc_nVa
        do bb = i_bb, f_bb
          if (bb == mb) cycle
          ba = bb - cc_nVa
          M2_A(ib,ja,aa,bb) = M2_A(ib,ja,aa,bb) &
          +1.0d0 * t2_A(ma, jb, ba, ab) * t2_A(ia, nb, na, mb) &
          +1.0d0 * t2_A(ia, nb, ba, ab) * t2_A(ma, jb, na, mb)
        enddo
      enddo
    enddo
  enddo

  ! ### Spin case: i_b, j_a, a_b, b_a ###

  do ib = i_ib, f_ib
    if (ib == nb) cycle
    ia = ib - cc_nOa
    do ja = i_ja, f_ja
      if (ja == ma) cycle
      jb = ja + cc_nOa
      do ab = i_ab, f_ab
        if (ab == mb) cycle
        aa = ab - cc_nVa
        do ba = i_ba, f_ba
          if (ba == na) cycle
          bb = ba + cc_nVa
          M2_A(ib,ja,ab,ba) = M2_A(ib,ja,ab,ba) &
          -1.0d0 * t2_A(ma, jb, aa, bb) * t2_A(ia, nb, na, mb) &
          -1.0d0 * t2_A(ia, nb, aa, bb) * t2_A(ma, jb, na, mb)
        enddo
      enddo
    enddo
  enddo

  ! ### Spin case: i_b, j_b, a_b, b_b ###

  do ib = i_ib, f_ib
    if (ib == nb) cycle
    ia = ib - cc_nOa
    do jb = i_jb, f_jb
      if (jb == nb) cycle
      ja = jb - cc_nOa
      do ab = i_ab, f_ab
        if (ab == mb) cycle
        aa = ab - cc_nVa
        do bb = i_bb, f_bb
          if (bb == mb) cycle
          ba = bb - cc_nVa
          M2_A(ib,jb,ab,bb) = M2_A(ib,jb,ab,bb) &
          -1.0d0 * t2_A(ma, ja, aa, ba) * t2_A(ia, nb, na, mb) &
          +1.0d0 * t2_A(ma, ia, aa, ba) * t2_A(ja, nb, na, mb)
        enddo
      enddo
    enddo
  enddo


end
