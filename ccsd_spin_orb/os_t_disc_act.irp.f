subroutine os_t_disc_act(nO,nV,det,t1_A,t2_A,M2_A)

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

  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      do ja = i_ja, f_ja
        if (ja == ma) cycle 
        jb = ja + cc_nOa
        do ia = i_ia, f_ia
          if (ia == ma) cycle 
          ib = ia + cc_nOa
          M2_A(ia,ja,aa,ba) = M2_A(ia,ja,aa,ba) & 
          +1.0d0 * t1_A(ia, aa) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) & 
          +1.0d0 * t1_A(ia, aa) * t1_A(ma, na) * t2_A(nb, jb, bb, mb) & 
          +1.0d0 * t1_A(ia, aa) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) & 
          -1.0d0 * t1_A(ia, aa) * t1_A(nb, mb) * t2_A(ma, jb, na, bb) & 
          +1.0d0 * t1_A(ia, aa) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(jb, mb) & 
          -1.0d0 * t1_A(ja, aa) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) & 
          -1.0d0 * t1_A(ja, aa) * t1_A(ma, na) * t2_A(nb, ib, bb, mb) & 
          -1.0d0 * t1_A(ja, aa) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) & 
          +1.0d0 * t1_A(ja, aa) * t1_A(nb, mb) * t2_A(ma, ib, na, bb) & 
          -1.0d0 * t1_A(ja, aa) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(ib, mb) & 
          -1.0d0 * t1_A(ia, ba) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) & 
          -1.0d0 * t1_A(ia, ba) * t1_A(ma, na) * t2_A(nb, jb, ab, mb) & 
          -1.0d0 * t1_A(ia, ba) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) & 
          +1.0d0 * t1_A(ia, ba) * t1_A(nb, mb) * t2_A(ma, jb, na, ab) & 
          -1.0d0 * t1_A(ia, ba) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(jb, mb) & 
          +1.0d0 * t1_A(ja, ba) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) & 
          +1.0d0 * t1_A(ja, ba) * t1_A(ma, na) * t2_A(nb, ib, ab, mb) & 
          +1.0d0 * t1_A(ja, ba) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) & 
          -1.0d0 * t1_A(ja, ba) * t1_A(nb, mb) * t2_A(ma, ib, na, ab) & 
          +1.0d0 * t1_A(ja, ba) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(ib, mb) & 
          -1.0d0 * t1_A(ib, ab) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) & 
          -1.0d0 * t1_A(ib, ab) * t1_A(ma, na) * t2_A(nb, jb, bb, mb) & 
          -1.0d0 * t1_A(ib, ab) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) & 
          +1.0d0 * t1_A(ib, ab) * t1_A(nb, mb) * t2_A(ma, jb, na, bb) & 
          -1.0d0 * t1_A(ib, ab) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(jb, mb) & 
          +1.0d0 * t1_A(jb, ab) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) & 
          +1.0d0 * t1_A(jb, ab) * t1_A(ma, na) * t2_A(nb, ib, bb, mb) & 
          +1.0d0 * t1_A(jb, ab) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) & 
          -1.0d0 * t1_A(jb, ab) * t1_A(nb, mb) * t2_A(ma, ib, na, bb) & 
          +1.0d0 * t1_A(jb, ab) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(ib, mb) & 
          +1.0d0 * t1_A(ib, bb) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) & 
          +1.0d0 * t1_A(ib, bb) * t1_A(ma, na) * t2_A(nb, jb, ab, mb) & 
          +1.0d0 * t1_A(ib, bb) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) & 
          -1.0d0 * t1_A(ib, bb) * t1_A(nb, mb) * t2_A(ma, jb, na, ab) & 
          +1.0d0 * t1_A(ib, bb) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(jb, mb) & 
          -1.0d0 * t1_A(jb, bb) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) & 
          -1.0d0 * t1_A(jb, bb) * t1_A(ma, na) * t2_A(nb, ib, ab, mb) & 
          -1.0d0 * t1_A(jb, bb) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) & 
          +1.0d0 * t1_A(jb, bb) * t1_A(nb, mb) * t2_A(ma, ib, na, ab) & 
          -1.0d0 * t1_A(jb, bb) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(ib, mb)
        enddo
      enddo
    enddo
  enddo

  !! Deltas:((na, ba))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    do ja = i_ja, f_ja
      if (ja == ma) cycle 
      jb = ja + cc_nOa
      do ia = i_ia, f_ia
        if (ia == ma) cycle 
        ib = ia + cc_nOa
        M2_A(ia,ja,aa,na) = M2_A(ia,ja,aa,na) & 
        +1.0d0 * t1_A(ia, aa) * t2_A(ma, jb, na, mb) & 
        +1.0d0 * t1_A(ia, aa) * t1_A(ma, na) * t1_A(jb, mb) & 
        -1.0d0 * t1_A(ja, aa) * t2_A(ma, ib, na, mb) & 
        -1.0d0 * t1_A(ja, aa) * t1_A(ma, na) * t1_A(ib, mb) & 
        -1.0d0 * t1_A(ia, na) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) & 
        -1.0d0 * t1_A(ia, na) * t1_A(ma, na) * t2_A(nb, jb, ab, mb) & 
        -1.0d0 * t1_A(ia, na) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) & 
        +1.0d0 * t1_A(ia, na) * t1_A(nb, mb) * t2_A(ma, jb, na, ab) & 
        -1.0d0 * t1_A(ia, na) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(jb, mb) & 
        +1.0d0 * t1_A(ja, na) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) & 
        +1.0d0 * t1_A(ja, na) * t1_A(ma, na) * t2_A(nb, ib, ab, mb) & 
        +1.0d0 * t1_A(ja, na) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) & 
        -1.0d0 * t1_A(ja, na) * t1_A(nb, mb) * t2_A(ma, ib, na, ab) & 
        +1.0d0 * t1_A(ja, na) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(ib, mb) & 
        -1.0d0 * t1_A(ib, ab) * t2_A(ma, jb, na, mb) & 
        -1.0d0 * t1_A(ib, ab) * t1_A(ma, na) * t1_A(jb, mb) & 
        +1.0d0 * t1_A(jb, ab) * t2_A(ma, ib, na, mb) & 
        +1.0d0 * t1_A(jb, ab) * t1_A(ma, na) * t1_A(ib, mb)
      enddo
    enddo
  enddo

  !! Deltas:((ma, ja))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      do ia = i_ia, f_ia
        if (ia == ma) cycle 
        ib = ia + cc_nOa
        M2_A(ia,ma,aa,ba) = M2_A(ia,ma,aa,ba) & 
        -1.0d0 * t1_A(ia, aa) * t2_A(ma, nb, na, bb) & 
        -1.0d0 * t1_A(ia, aa) * t1_A(nb, bb) * t1_A(ma, na) & 
        -1.0d0 * t1_A(ma, aa) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) & 
        -1.0d0 * t1_A(ma, aa) * t1_A(ma, na) * t2_A(nb, ib, bb, mb) & 
        -1.0d0 * t1_A(ma, aa) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) & 
        +1.0d0 * t1_A(ma, aa) * t1_A(nb, mb) * t2_A(ma, ib, na, bb) & 
        -1.0d0 * t1_A(ma, aa) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(ib, mb) & 
        +1.0d0 * t1_A(ia, ba) * t2_A(ma, nb, na, ab) & 
        +1.0d0 * t1_A(ia, ba) * t1_A(nb, ab) * t1_A(ma, na) & 
        +1.0d0 * t1_A(ma, ba) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) & 
        +1.0d0 * t1_A(ma, ba) * t1_A(ma, na) * t2_A(nb, ib, ab, mb) & 
        +1.0d0 * t1_A(ma, ba) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) & 
        -1.0d0 * t1_A(ma, ba) * t1_A(nb, mb) * t2_A(ma, ib, na, ab) & 
        +1.0d0 * t1_A(ma, ba) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(ib, mb) & 
        +1.0d0 * t1_A(ib, ab) * t2_A(ma, nb, na, bb) & 
        +1.0d0 * t1_A(ib, ab) * t1_A(nb, bb) * t1_A(ma, na) & 
        -1.0d0 * t1_A(ib, bb) * t2_A(ma, nb, na, ab) & 
        -1.0d0 * t1_A(ib, bb) * t1_A(nb, ab) * t1_A(ma, na)
      enddo
    enddo
  enddo

  !! Deltas:((na, aa))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ja = i_ja, f_ja
      if (ja == ma) cycle 
      jb = ja + cc_nOa
      do ia = i_ia, f_ia
        if (ia == ma) cycle 
        ib = ia + cc_nOa
        M2_A(ia,ja,na,ba) = M2_A(ia,ja,na,ba) & 
        +1.0d0 * t1_A(ia, na) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) & 
        +1.0d0 * t1_A(ia, na) * t1_A(ma, na) * t2_A(nb, jb, bb, mb) & 
        +1.0d0 * t1_A(ia, na) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) & 
        -1.0d0 * t1_A(ia, na) * t1_A(nb, mb) * t2_A(ma, jb, na, bb) & 
        +1.0d0 * t1_A(ia, na) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(jb, mb) & 
        -1.0d0 * t1_A(ja, na) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) & 
        -1.0d0 * t1_A(ja, na) * t1_A(ma, na) * t2_A(nb, ib, bb, mb) & 
        -1.0d0 * t1_A(ja, na) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) & 
        +1.0d0 * t1_A(ja, na) * t1_A(nb, mb) * t2_A(ma, ib, na, bb) & 
        -1.0d0 * t1_A(ja, na) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(ib, mb) & 
        -1.0d0 * t1_A(ia, ba) * t2_A(ma, jb, na, mb) & 
        -1.0d0 * t1_A(ia, ba) * t1_A(ma, na) * t1_A(jb, mb) & 
        +1.0d0 * t1_A(ja, ba) * t2_A(ma, ib, na, mb) & 
        +1.0d0 * t1_A(ja, ba) * t1_A(ma, na) * t1_A(ib, mb) & 
        +1.0d0 * t1_A(ib, bb) * t2_A(ma, jb, na, mb) & 
        +1.0d0 * t1_A(ib, bb) * t1_A(ma, na) * t1_A(jb, mb) & 
        -1.0d0 * t1_A(jb, bb) * t2_A(ma, ib, na, mb) & 
        -1.0d0 * t1_A(jb, bb) * t1_A(ma, na) * t1_A(ib, mb)
      enddo
    enddo
  enddo

  !! Deltas:((ma, ia))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      do ja = i_ja, f_ja
        if (ja == ma) cycle 
        jb = ja + cc_nOa
        M2_A(ma,ja,aa,ba) = M2_A(ma,ja,aa,ba) & 
        +1.0d0 * t1_A(ma, aa) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) & 
        +1.0d0 * t1_A(ma, aa) * t1_A(ma, na) * t2_A(nb, jb, bb, mb) & 
        +1.0d0 * t1_A(ma, aa) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) & 
        -1.0d0 * t1_A(ma, aa) * t1_A(nb, mb) * t2_A(ma, jb, na, bb) & 
        +1.0d0 * t1_A(ma, aa) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(jb, mb) & 
        +1.0d0 * t1_A(ja, aa) * t2_A(ma, nb, na, bb) & 
        +1.0d0 * t1_A(ja, aa) * t1_A(nb, bb) * t1_A(ma, na) & 
        -1.0d0 * t1_A(ma, ba) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) & 
        -1.0d0 * t1_A(ma, ba) * t1_A(ma, na) * t2_A(nb, jb, ab, mb) & 
        -1.0d0 * t1_A(ma, ba) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) & 
        +1.0d0 * t1_A(ma, ba) * t1_A(nb, mb) * t2_A(ma, jb, na, ab) & 
        -1.0d0 * t1_A(ma, ba) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(jb, mb) & 
        -1.0d0 * t1_A(ja, ba) * t2_A(ma, nb, na, ab) & 
        -1.0d0 * t1_A(ja, ba) * t1_A(nb, ab) * t1_A(ma, na) & 
        -1.0d0 * t1_A(jb, ab) * t2_A(ma, nb, na, bb) & 
        -1.0d0 * t1_A(jb, ab) * t1_A(nb, bb) * t1_A(ma, na) & 
        +1.0d0 * t1_A(jb, bb) * t2_A(ma, nb, na, ab) & 
        +1.0d0 * t1_A(jb, bb) * t1_A(nb, ab) * t1_A(ma, na)
      enddo
    enddo
  enddo

  !! Deltas:((ma, ja), (na, ba))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    do ia = i_ia, f_ia
      if (ia == ma) cycle 
      ib = ia + cc_nOa
      M2_A(ia,ma,aa,na) = M2_A(ia,ma,aa,na) & 
      -1.0d0 * t1_A(ia, aa) * t1_A(ma, na) & 
      -1.0d0 * t1_A(ma, aa) * t2_A(ma, ib, na, mb) & 
      -1.0d0 * t1_A(ma, aa) * t1_A(ma, na) * t1_A(ib, mb) & 
      +1.0d0 * t1_A(ia, na) * t2_A(ma, nb, na, ab) & 
      +1.0d0 * t1_A(ia, na) * t1_A(nb, ab) * t1_A(ma, na) & 
      +1.0d0 * t1_A(ma, na) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(ma, na) * t2_A(nb, ib, ab, mb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) & 
      -1.0d0 * t1_A(ma, na) * t1_A(nb, mb) * t2_A(ma, ib, na, ab) & 
      +1.0d0 * t1_A(ma, na) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(ib, mb) & 
      +1.0d0 * t1_A(ib, ab) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((ma, ja), (na, aa))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ia = i_ia, f_ia
      if (ia == ma) cycle 
      ib = ia + cc_nOa
      M2_A(ia,ma,na,ba) = M2_A(ia,ma,na,ba) & 
      -1.0d0 * t1_A(ia, na) * t2_A(ma, nb, na, bb) & 
      -1.0d0 * t1_A(ia, na) * t1_A(nb, bb) * t1_A(ma, na) & 
      -1.0d0 * t1_A(ma, na) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(ma, na) * t2_A(nb, ib, bb, mb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(nb, mb) * t2_A(ma, ib, na, bb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(ib, mb) & 
      +1.0d0 * t1_A(ia, ba) * t1_A(ma, na) & 
      +1.0d0 * t1_A(ma, ba) * t2_A(ma, ib, na, mb) & 
      +1.0d0 * t1_A(ma, ba) * t1_A(ma, na) * t1_A(ib, mb) & 
      -1.0d0 * t1_A(ib, bb) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((ma, ia), (na, ba))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    do ja = i_ja, f_ja
      if (ja == ma) cycle 
      jb = ja + cc_nOa
      M2_A(ma,ja,aa,na) = M2_A(ma,ja,aa,na) & 
      +1.0d0 * t1_A(ma, aa) * t2_A(ma, jb, na, mb) & 
      +1.0d0 * t1_A(ma, aa) * t1_A(ma, na) * t1_A(jb, mb) & 
      +1.0d0 * t1_A(ja, aa) * t1_A(ma, na) & 
      -1.0d0 * t1_A(ma, na) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(ma, na) * t2_A(nb, jb, ab, mb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) & 
      +1.0d0 * t1_A(ma, na) * t1_A(nb, mb) * t2_A(ma, jb, na, ab) & 
      -1.0d0 * t1_A(ma, na) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(jb, mb) & 
      -1.0d0 * t1_A(ja, na) * t2_A(ma, nb, na, ab) & 
      -1.0d0 * t1_A(ja, na) * t1_A(nb, ab) * t1_A(ma, na) & 
      -1.0d0 * t1_A(jb, ab) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((ma, ia), (na, aa))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ja = i_ja, f_ja
      if (ja == ma) cycle 
      jb = ja + cc_nOa
      M2_A(ma,ja,na,ba) = M2_A(ma,ja,na,ba) & 
      +1.0d0 * t1_A(ma, na) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(ma, na) * t2_A(nb, jb, bb, mb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(nb, mb) * t2_A(ma, jb, na, bb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(jb, mb) & 
      +1.0d0 * t1_A(ja, na) * t2_A(ma, nb, na, bb) & 
      +1.0d0 * t1_A(ja, na) * t1_A(nb, bb) * t1_A(ma, na) & 
      -1.0d0 * t1_A(ma, ba) * t2_A(ma, jb, na, mb) & 
      -1.0d0 * t1_A(ma, ba) * t1_A(ma, na) * t1_A(jb, mb) & 
      -1.0d0 * t1_A(ja, ba) * t1_A(ma, na) & 
      +1.0d0 * t1_A(jb, bb) * t1_A(ma, na)
    enddo
  enddo

  ! ### Spin case: i_a, j_b, a_a, b_b ###

  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      do jb = i_jb, f_jb
        if (jb == nb) cycle 
        ja = jb - cc_nOa
        do ia = i_ia, f_ia
          if (ia == ma) cycle 
          ib = ia + cc_nOa
          M2_A(ia,jb,aa,bb) = M2_A(ia,jb,aa,bb) & 
          +1.0d0 * t1_A(ia, aa) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) & 
          +1.0d0 * t1_A(ia, aa) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) & 
          -1.0d0 * t1_A(ia, aa) * t1_A(ma, na) * t2_A(ja, nb, ba, mb) & 
          +1.0d0 * t1_A(ia, aa) * t1_A(nb, mb) * t2_A(ma, ja, ba, na) & 
          +1.0d0 * t1_A(ia, aa) * t1_A(ma, ba) * t1_A(ja, na) * t1_A(nb, mb) & 
          +1.0d0 * t1_A(jb, bb) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) & 
          +1.0d0 * t1_A(jb, bb) * t1_A(ma, na) * t2_A(nb, ib, ab, mb) & 
          +1.0d0 * t1_A(jb, bb) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) & 
          -1.0d0 * t1_A(jb, bb) * t1_A(nb, mb) * t2_A(ma, ib, na, ab) & 
          +1.0d0 * t1_A(jb, bb) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(ib, mb) & 
          -1.0d0 * t1_A(ib, ab) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) & 
          -1.0d0 * t1_A(ib, ab) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) & 
          +1.0d0 * t1_A(ib, ab) * t1_A(ma, na) * t2_A(ja, nb, ba, mb) & 
          -1.0d0 * t1_A(ib, ab) * t1_A(nb, mb) * t2_A(ma, ja, ba, na) & 
          -1.0d0 * t1_A(ib, ab) * t1_A(ma, ba) * t1_A(ja, na) * t1_A(nb, mb) & 
          -1.0d0 * t1_A(ja, ba) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) & 
          -1.0d0 * t1_A(ja, ba) * t1_A(ma, na) * t2_A(nb, ib, ab, mb) & 
          -1.0d0 * t1_A(ja, ba) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) & 
          +1.0d0 * t1_A(ja, ba) * t1_A(nb, mb) * t2_A(ma, ib, na, ab) & 
          -1.0d0 * t1_A(ja, ba) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(ib, mb)
        enddo
      enddo
    enddo
  enddo

  !! Deltas:((mb, bb))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    do jb = i_jb, f_jb
      if (jb == nb) cycle 
      ja = jb - cc_nOa
      do ia = i_ia, f_ia
        if (ia == ma) cycle 
        ib = ia + cc_nOa
        M2_A(ia,jb,aa,mb) = M2_A(ia,jb,aa,mb) & 
        +1.0d0 * t1_A(ia, aa) * t2_A(ja, nb, na, mb) & 
        +1.0d0 * t1_A(ia, aa) * t1_A(ja, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(jb, mb) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) & 
        +1.0d0 * t1_A(jb, mb) * t1_A(ma, na) * t2_A(nb, ib, ab, mb) & 
        +1.0d0 * t1_A(jb, mb) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) & 
        -1.0d0 * t1_A(jb, mb) * t1_A(nb, mb) * t2_A(ma, ib, na, ab) & 
        +1.0d0 * t1_A(jb, mb) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(ib, mb) & 
        -1.0d0 * t1_A(ib, ab) * t2_A(ja, nb, na, mb) & 
        -1.0d0 * t1_A(ib, ab) * t1_A(ja, na) * t1_A(nb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((nb, jb))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      do ia = i_ia, f_ia
        if (ia == ma) cycle 
        ib = ia + cc_nOa
        M2_A(ia,nb,aa,bb) = M2_A(ia,nb,aa,bb) & 
        -1.0d0 * t1_A(ia, aa) * t2_A(ma, nb, ba, mb) & 
        -1.0d0 * t1_A(ia, aa) * t1_A(ma, ba) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(nb, bb) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) & 
        +1.0d0 * t1_A(nb, bb) * t1_A(ma, na) * t2_A(nb, ib, ab, mb) & 
        +1.0d0 * t1_A(nb, bb) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) & 
        -1.0d0 * t1_A(nb, bb) * t1_A(nb, mb) * t2_A(ma, ib, na, ab) & 
        +1.0d0 * t1_A(nb, bb) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(ib, mb) & 
        +1.0d0 * t1_A(ib, ab) * t2_A(ma, nb, ba, mb) & 
        +1.0d0 * t1_A(ib, ab) * t1_A(ma, ba) * t1_A(nb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((na, aa))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do jb = i_jb, f_jb
      if (jb == nb) cycle 
      ja = jb - cc_nOa
      do ia = i_ia, f_ia
        if (ia == ma) cycle 
        ib = ia + cc_nOa
        M2_A(ia,jb,na,bb) = M2_A(ia,jb,na,bb) & 
        +1.0d0 * t1_A(ia, na) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) & 
        +1.0d0 * t1_A(ia, na) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) & 
        -1.0d0 * t1_A(ia, na) * t1_A(ma, na) * t2_A(ja, nb, ba, mb) & 
        +1.0d0 * t1_A(ia, na) * t1_A(nb, mb) * t2_A(ma, ja, ba, na) & 
        +1.0d0 * t1_A(ia, na) * t1_A(ma, ba) * t1_A(ja, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(jb, bb) * t2_A(ma, ib, na, mb) & 
        +1.0d0 * t1_A(jb, bb) * t1_A(ma, na) * t1_A(ib, mb) & 
        -1.0d0 * t1_A(ja, ba) * t2_A(ma, ib, na, mb) & 
        -1.0d0 * t1_A(ja, ba) * t1_A(ma, na) * t1_A(ib, mb)
      enddo
    enddo
  enddo

  !! Deltas:((ma, ia))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      do jb = i_jb, f_jb
        if (jb == nb) cycle 
        ja = jb - cc_nOa
        M2_A(ma,jb,aa,bb) = M2_A(ma,jb,aa,bb) & 
        +1.0d0 * t1_A(ma, aa) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) & 
        +1.0d0 * t1_A(ma, aa) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) & 
        -1.0d0 * t1_A(ma, aa) * t1_A(ma, na) * t2_A(ja, nb, ba, mb) & 
        +1.0d0 * t1_A(ma, aa) * t1_A(nb, mb) * t2_A(ma, ja, ba, na) & 
        +1.0d0 * t1_A(ma, aa) * t1_A(ma, ba) * t1_A(ja, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(jb, bb) * t2_A(ma, nb, na, ab) & 
        -1.0d0 * t1_A(jb, bb) * t1_A(nb, ab) * t1_A(ma, na) & 
        +1.0d0 * t1_A(ja, ba) * t2_A(ma, nb, na, ab) & 
        +1.0d0 * t1_A(ja, ba) * t1_A(nb, ab) * t1_A(ma, na)
      enddo
    enddo
  enddo

  !! Deltas:((nb, jb), (mb, bb))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    do ia = i_ia, f_ia
      if (ia == ma) cycle 
      ib = ia + cc_nOa
      M2_A(ia,nb,aa,mb) = M2_A(ia,nb,aa,mb) & 
      -1.0d0 * t1_A(ia, aa) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(ma, na) * t2_A(nb, ib, ab, mb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(nb, mb) * t2_A(ma, ib, na, ab) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(ib, mb) & 
      +1.0d0 * t1_A(ib, ab) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((mb, bb), (na, aa))
  do jb = i_jb, f_jb
    if (jb == nb) cycle 
    ja = jb - cc_nOa
    do ia = i_ia, f_ia
      if (ia == ma) cycle 
      ib = ia + cc_nOa
      M2_A(ia,jb,na,mb) = M2_A(ia,jb,na,mb) & 
      +1.0d0 * t1_A(ia, na) * t2_A(ja, nb, na, mb) & 
      +1.0d0 * t1_A(ia, na) * t1_A(ja, na) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((nb, jb), (na, aa))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do ia = i_ia, f_ia
      if (ia == ma) cycle 
      ib = ia + cc_nOa
      M2_A(ia,nb,na,bb) = M2_A(ia,nb,na,bb) & 
      -1.0d0 * t1_A(ia, na) * t2_A(ma, nb, ba, mb) & 
      -1.0d0 * t1_A(ia, na) * t1_A(ma, ba) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(nb, bb) * t2_A(ma, ib, na, mb) & 
      +1.0d0 * t1_A(nb, bb) * t1_A(ma, na) * t1_A(ib, mb)
    enddo
  enddo

  !! Deltas:((ma, ia), (mb, bb))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    do jb = i_jb, f_jb
      if (jb == nb) cycle 
      ja = jb - cc_nOa
      M2_A(ma,jb,aa,mb) = M2_A(ma,jb,aa,mb) & 
      +1.0d0 * t1_A(ma, aa) * t2_A(ja, nb, na, mb) & 
      +1.0d0 * t1_A(ma, aa) * t1_A(ja, na) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(jb, mb) * t2_A(ma, nb, na, ab) & 
      -1.0d0 * t1_A(jb, mb) * t1_A(nb, ab) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((ma, ia), (nb, jb))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      M2_A(ma,nb,aa,bb) = M2_A(ma,nb,aa,bb) & 
      -1.0d0 * t1_A(ma, aa) * t2_A(ma, nb, ba, mb) & 
      -1.0d0 * t1_A(ma, aa) * t1_A(ma, ba) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((ma, ia), (na, aa))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do jb = i_jb, f_jb
      if (jb == nb) cycle 
      ja = jb - cc_nOa
      M2_A(ma,jb,na,bb) = M2_A(ma,jb,na,bb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(ma, na) * t2_A(ja, nb, ba, mb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(nb, mb) * t2_A(ma, ja, ba, na) & 
      +1.0d0 * t1_A(ma, na) * t1_A(ma, ba) * t1_A(ja, na) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(jb, bb) * t1_A(ma, na) & 
      +1.0d0 * t1_A(ja, ba) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((na, aa), (mb, bb))
  do jb = i_jb, f_jb
    if (jb == nb) cycle 
    ja = jb - cc_nOa
    do ia = i_ia, f_ia
      if (ia == ma) cycle 
      ib = ia + cc_nOa
      M2_A(ia,jb,na,mb) = M2_A(ia,jb,na,mb) & 
      +1.0d0 * t1_A(jb, mb) * t2_A(ma, ib, na, mb) & 
      +1.0d0 * t1_A(jb, mb) * t1_A(ma, na) * t1_A(ib, mb)
    enddo
  enddo

  !! Deltas:((nb, jb), (ma, ia))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      M2_A(ma,nb,aa,bb) = M2_A(ma,nb,aa,bb) & 
      -1.0d0 * t1_A(nb, bb) * t2_A(ma, nb, na, ab) & 
      -1.0d0 * t1_A(nb, bb) * t1_A(nb, ab) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((nb, jb), (na, aa), (mb, bb))
  do ia = i_ia, f_ia
    if (ia == ma) cycle 
    ib = ia + cc_nOa
    M2_A(ia,nb,na,mb) = M2_A(ia,nb,na,mb) & 
    -1.0d0 * t1_A(ia, na) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ia), (mb, bb), (nb, jb))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    M2_A(ma,nb,aa,mb) = M2_A(ma,nb,aa,mb) & 
    -1.0d0 * t1_A(ma, aa) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ia), (na, aa), (mb, bb))
  do jb = i_jb, f_jb
    if (jb == nb) cycle 
    ja = jb - cc_nOa
    M2_A(ma,jb,na,mb) = M2_A(ma,jb,na,mb) & 
    +1.0d0 * t1_A(ma, na) * t2_A(ja, nb, na, mb) & 
    +1.0d0 * t1_A(ma, na) * t1_A(ja, na) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ia), (na, aa), (nb, jb))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    M2_A(ma,nb,na,bb) = M2_A(ma,nb,na,bb) & 
    -1.0d0 * t1_A(ma, na) * t2_A(ma, nb, ba, mb) & 
    -1.0d0 * t1_A(ma, na) * t1_A(ma, ba) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ia), (mb, bb), (na, aa))
  do jb = i_jb, f_jb
    if (jb == nb) cycle 
    ja = jb - cc_nOa
    M2_A(ma,jb,na,mb) = M2_A(ma,jb,na,mb) & 
    -1.0d0 * t1_A(jb, mb) * t1_A(ma, na)
  enddo

  !! Deltas:((nb, jb), (na, aa), (ma, ia))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    M2_A(ma,nb,na,bb) = M2_A(ma,nb,na,bb) & 
    -1.0d0 * t1_A(nb, bb) * t1_A(ma, na)
  enddo

  !! Deltas:((nb, jb), (mb, bb), (na, aa))
  do ia = i_ia, f_ia
    if (ia == ma) cycle 
    ib = ia + cc_nOa
    M2_A(ia,nb,na,mb) = M2_A(ia,nb,na,mb) & 
    +1.0d0 * t1_A(nb, mb) * t2_A(ma, ib, na, mb) & 
    +1.0d0 * t1_A(nb, mb) * t1_A(ma, na) * t1_A(ib, mb)
  enddo

  !! Deltas:((nb, jb), (mb, bb), (ma, ia))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    M2_A(ma,nb,aa,mb) = M2_A(ma,nb,aa,mb) & 
    -1.0d0 * t1_A(nb, mb) * t2_A(ma, nb, na, ab) & 
    -1.0d0 * t1_A(nb, mb) * t1_A(nb, ab) * t1_A(ma, na)
  enddo

  !! Deltas:((ma, ia), (na, aa), (mb, bb), (nb, jb))
  M2_A(ma,nb,na,mb) = M2_A(ma,nb,na,mb) & 
  -1.0d0 * t1_A(ma, na) * t1_A(nb, mb)

  !! Deltas:((nb, jb), (mb, bb), (na, aa), (ma, ia))
  M2_A(ma,nb,na,mb) = M2_A(ma,nb,na,mb) & 
  -1.0d0 * t1_A(nb, mb) * t1_A(ma, na)

  ! ### Spin case: i_a, j_b, a_b, b_a ###

  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      do jb = i_jb, f_jb
        if (jb == nb) cycle 
        ja = jb - cc_nOa
        do ia = i_ia, f_ia
          if (ia == ma) cycle 
          ib = ia + cc_nOa
          M2_A(ia,jb,ab,ba) = M2_A(ia,jb,ab,ba) & 
          -1.0d0 * t1_A(jb, ab) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) & 
          -1.0d0 * t1_A(jb, ab) * t1_A(ma, na) * t2_A(nb, ib, bb, mb) & 
          -1.0d0 * t1_A(jb, ab) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) & 
          +1.0d0 * t1_A(jb, ab) * t1_A(nb, mb) * t2_A(ma, ib, na, bb) & 
          -1.0d0 * t1_A(jb, ab) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(ib, mb) & 
          -1.0d0 * t1_A(ia, ba) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) & 
          -1.0d0 * t1_A(ia, ba) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) & 
          +1.0d0 * t1_A(ia, ba) * t1_A(ma, na) * t2_A(ja, nb, aa, mb) & 
          -1.0d0 * t1_A(ia, ba) * t1_A(nb, mb) * t2_A(ma, ja, aa, na) & 
          -1.0d0 * t1_A(ia, ba) * t1_A(ma, aa) * t1_A(ja, na) * t1_A(nb, mb) & 
          +1.0d0 * t1_A(ja, aa) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) & 
          +1.0d0 * t1_A(ja, aa) * t1_A(ma, na) * t2_A(nb, ib, bb, mb) & 
          +1.0d0 * t1_A(ja, aa) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) & 
          -1.0d0 * t1_A(ja, aa) * t1_A(nb, mb) * t2_A(ma, ib, na, bb) & 
          +1.0d0 * t1_A(ja, aa) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(ib, mb) & 
          +1.0d0 * t1_A(ib, bb) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) & 
          +1.0d0 * t1_A(ib, bb) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) & 
          -1.0d0 * t1_A(ib, bb) * t1_A(ma, na) * t2_A(ja, nb, aa, mb) & 
          +1.0d0 * t1_A(ib, bb) * t1_A(nb, mb) * t2_A(ma, ja, aa, na) & 
          +1.0d0 * t1_A(ib, bb) * t1_A(ma, aa) * t1_A(ja, na) * t1_A(nb, mb)
        enddo
      enddo
    enddo
  enddo

  !! Deltas:((na, ba))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    do jb = i_jb, f_jb
      if (jb == nb) cycle 
      ja = jb - cc_nOa
      do ia = i_ia, f_ia
        if (ia == ma) cycle 
        ib = ia + cc_nOa
        M2_A(ia,jb,ab,na) = M2_A(ia,jb,ab,na) & 
        -1.0d0 * t1_A(jb, ab) * t2_A(ma, ib, na, mb) & 
        -1.0d0 * t1_A(jb, ab) * t1_A(ma, na) * t1_A(ib, mb) & 
        -1.0d0 * t1_A(ia, na) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) & 
        -1.0d0 * t1_A(ia, na) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) & 
        +1.0d0 * t1_A(ia, na) * t1_A(ma, na) * t2_A(ja, nb, aa, mb) & 
        -1.0d0 * t1_A(ia, na) * t1_A(nb, mb) * t2_A(ma, ja, aa, na) & 
        -1.0d0 * t1_A(ia, na) * t1_A(ma, aa) * t1_A(ja, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(ja, aa) * t2_A(ma, ib, na, mb) & 
        +1.0d0 * t1_A(ja, aa) * t1_A(ma, na) * t1_A(ib, mb)
      enddo
    enddo
  enddo

  !! Deltas:((ma, ia))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      do jb = i_jb, f_jb
        if (jb == nb) cycle 
        ja = jb - cc_nOa
        M2_A(ma,jb,ab,ba) = M2_A(ma,jb,ab,ba) & 
        +1.0d0 * t1_A(jb, ab) * t2_A(ma, nb, na, bb) & 
        +1.0d0 * t1_A(jb, ab) * t1_A(nb, bb) * t1_A(ma, na) & 
        -1.0d0 * t1_A(ma, ba) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) & 
        -1.0d0 * t1_A(ma, ba) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) & 
        +1.0d0 * t1_A(ma, ba) * t1_A(ma, na) * t2_A(ja, nb, aa, mb) & 
        -1.0d0 * t1_A(ma, ba) * t1_A(nb, mb) * t2_A(ma, ja, aa, na) & 
        -1.0d0 * t1_A(ma, ba) * t1_A(ma, aa) * t1_A(ja, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(ja, aa) * t2_A(ma, nb, na, bb) & 
        -1.0d0 * t1_A(ja, aa) * t1_A(nb, bb) * t1_A(ma, na)
      enddo
    enddo
  enddo

  !! Deltas:((mb, ab))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do jb = i_jb, f_jb
      if (jb == nb) cycle 
      ja = jb - cc_nOa
      do ia = i_ia, f_ia
        if (ia == ma) cycle 
        ib = ia + cc_nOa
        M2_A(ia,jb,mb,ba) = M2_A(ia,jb,mb,ba) & 
        -1.0d0 * t1_A(jb, mb) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) & 
        -1.0d0 * t1_A(jb, mb) * t1_A(ma, na) * t2_A(nb, ib, bb, mb) & 
        -1.0d0 * t1_A(jb, mb) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) & 
        +1.0d0 * t1_A(jb, mb) * t1_A(nb, mb) * t2_A(ma, ib, na, bb) & 
        -1.0d0 * t1_A(jb, mb) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(ib, mb) & 
        -1.0d0 * t1_A(ia, ba) * t2_A(ja, nb, na, mb) & 
        -1.0d0 * t1_A(ia, ba) * t1_A(ja, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(ib, bb) * t2_A(ja, nb, na, mb) & 
        +1.0d0 * t1_A(ib, bb) * t1_A(ja, na) * t1_A(nb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((nb, jb))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      do ia = i_ia, f_ia
        if (ia == ma) cycle 
        ib = ia + cc_nOa
        M2_A(ia,nb,ab,ba) = M2_A(ia,nb,ab,ba) & 
        -1.0d0 * t1_A(nb, ab) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) & 
        -1.0d0 * t1_A(nb, ab) * t1_A(ma, na) * t2_A(nb, ib, bb, mb) & 
        -1.0d0 * t1_A(nb, ab) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) & 
        +1.0d0 * t1_A(nb, ab) * t1_A(nb, mb) * t2_A(ma, ib, na, bb) & 
        -1.0d0 * t1_A(nb, ab) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(ib, mb) & 
        +1.0d0 * t1_A(ia, ba) * t2_A(ma, nb, aa, mb) & 
        +1.0d0 * t1_A(ia, ba) * t1_A(ma, aa) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(ib, bb) * t2_A(ma, nb, aa, mb) & 
        -1.0d0 * t1_A(ib, bb) * t1_A(ma, aa) * t1_A(nb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((ma, ia), (na, ba))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    do jb = i_jb, f_jb
      if (jb == nb) cycle 
      ja = jb - cc_nOa
      M2_A(ma,jb,ab,na) = M2_A(ma,jb,ab,na) & 
      +1.0d0 * t1_A(jb, ab) * t1_A(ma, na) & 
      -1.0d0 * t1_A(ma, na) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(ma, na) * t2_A(ja, nb, aa, mb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(nb, mb) * t2_A(ma, ja, aa, na) & 
      -1.0d0 * t1_A(ma, na) * t1_A(ma, aa) * t1_A(ja, na) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(ja, aa) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((na, ba), (mb, ab))
  do jb = i_jb, f_jb
    if (jb == nb) cycle 
    ja = jb - cc_nOa
    do ia = i_ia, f_ia
      if (ia == ma) cycle 
      ib = ia + cc_nOa
      M2_A(ia,jb,mb,na) = M2_A(ia,jb,mb,na) & 
      -1.0d0 * t1_A(jb, mb) * t2_A(ma, ib, na, mb) & 
      -1.0d0 * t1_A(jb, mb) * t1_A(ma, na) * t1_A(ib, mb)
    enddo
  enddo

  !! Deltas:((ma, ia), (mb, ab))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do jb = i_jb, f_jb
      if (jb == nb) cycle 
      ja = jb - cc_nOa
      M2_A(ma,jb,mb,ba) = M2_A(ma,jb,mb,ba) & 
      +1.0d0 * t1_A(jb, mb) * t2_A(ma, nb, na, bb) & 
      +1.0d0 * t1_A(jb, mb) * t1_A(nb, bb) * t1_A(ma, na) & 
      -1.0d0 * t1_A(ma, ba) * t2_A(ja, nb, na, mb) & 
      -1.0d0 * t1_A(ma, ba) * t1_A(ja, na) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((nb, jb), (na, ba))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    do ia = i_ia, f_ia
      if (ia == ma) cycle 
      ib = ia + cc_nOa
      M2_A(ia,nb,ab,na) = M2_A(ia,nb,ab,na) & 
      -1.0d0 * t1_A(nb, ab) * t2_A(ma, ib, na, mb) & 
      -1.0d0 * t1_A(nb, ab) * t1_A(ma, na) * t1_A(ib, mb) & 
      +1.0d0 * t1_A(ia, na) * t2_A(ma, nb, aa, mb) & 
      +1.0d0 * t1_A(ia, na) * t1_A(ma, aa) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((nb, jb), (ma, ia))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      M2_A(ma,nb,ab,ba) = M2_A(ma,nb,ab,ba) & 
      +1.0d0 * t1_A(nb, ab) * t2_A(ma, nb, na, bb) & 
      +1.0d0 * t1_A(nb, ab) * t1_A(nb, bb) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((nb, jb), (mb, ab))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ia = i_ia, f_ia
      if (ia == ma) cycle 
      ib = ia + cc_nOa
      M2_A(ia,nb,mb,ba) = M2_A(ia,nb,mb,ba) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(ma, na) * t2_A(nb, ib, bb, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(nb, mb) * t2_A(ma, ib, na, bb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(ib, mb) & 
      +1.0d0 * t1_A(ia, ba) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(ib, bb) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((mb, ab), (na, ba))
  do jb = i_jb, f_jb
    if (jb == nb) cycle 
    ja = jb - cc_nOa
    do ia = i_ia, f_ia
      if (ia == ma) cycle 
      ib = ia + cc_nOa
      M2_A(ia,jb,mb,na) = M2_A(ia,jb,mb,na) & 
      -1.0d0 * t1_A(ia, na) * t2_A(ja, nb, na, mb) & 
      -1.0d0 * t1_A(ia, na) * t1_A(ja, na) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((ma, ia), (nb, jb))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      M2_A(ma,nb,ab,ba) = M2_A(ma,nb,ab,ba) & 
      +1.0d0 * t1_A(ma, ba) * t2_A(ma, nb, aa, mb) & 
      +1.0d0 * t1_A(ma, ba) * t1_A(ma, aa) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((ma, ia), (mb, ab), (na, ba))
  do jb = i_jb, f_jb
    if (jb == nb) cycle 
    ja = jb - cc_nOa
    M2_A(ma,jb,mb,na) = M2_A(ma,jb,mb,na) & 
    +1.0d0 * t1_A(jb, mb) * t1_A(ma, na)
  enddo

  !! Deltas:((nb, jb), (na, ba), (ma, ia))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    M2_A(ma,nb,ab,na) = M2_A(ma,nb,ab,na) & 
    +1.0d0 * t1_A(nb, ab) * t1_A(ma, na)
  enddo

  !! Deltas:((nb, jb), (mb, ab), (na, ba))
  do ia = i_ia, f_ia
    if (ia == ma) cycle 
    ib = ia + cc_nOa
    M2_A(ia,nb,mb,na) = M2_A(ia,nb,mb,na) & 
    -1.0d0 * t1_A(nb, mb) * t2_A(ma, ib, na, mb) & 
    -1.0d0 * t1_A(nb, mb) * t1_A(ma, na) * t1_A(ib, mb)
  enddo

  !! Deltas:((nb, jb), (mb, ab), (ma, ia))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    M2_A(ma,nb,mb,ba) = M2_A(ma,nb,mb,ba) & 
    +1.0d0 * t1_A(nb, mb) * t2_A(ma, nb, na, bb) & 
    +1.0d0 * t1_A(nb, mb) * t1_A(nb, bb) * t1_A(ma, na)
  enddo

  !! Deltas:((nb, jb), (na, ba), (mb, ab))
  do ia = i_ia, f_ia
    if (ia == ma) cycle 
    ib = ia + cc_nOa
    M2_A(ia,nb,mb,na) = M2_A(ia,nb,mb,na) & 
    +1.0d0 * t1_A(ia, na) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ia), (mb, ab), (nb, jb))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    M2_A(ma,nb,mb,ba) = M2_A(ma,nb,mb,ba) & 
    +1.0d0 * t1_A(ma, ba) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ia), (na, ba), (mb, ab))
  do jb = i_jb, f_jb
    if (jb == nb) cycle 
    ja = jb - cc_nOa
    M2_A(ma,jb,mb,na) = M2_A(ma,jb,mb,na) & 
    -1.0d0 * t1_A(ma, na) * t2_A(ja, nb, na, mb) & 
    -1.0d0 * t1_A(ma, na) * t1_A(ja, na) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ia), (na, ba), (nb, jb))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    M2_A(ma,nb,ab,na) = M2_A(ma,nb,ab,na) & 
    +1.0d0 * t1_A(ma, na) * t2_A(ma, nb, aa, mb) & 
    +1.0d0 * t1_A(ma, na) * t1_A(ma, aa) * t1_A(nb, mb)
  enddo

  !! Deltas:((nb, jb), (mb, ab), (na, ba), (ma, ia))
  M2_A(ma,nb,mb,na) = M2_A(ma,nb,mb,na) & 
  +1.0d0 * t1_A(nb, mb) * t1_A(ma, na)

  !! Deltas:((ma, ia), (na, ba), (mb, ab), (nb, jb))
  M2_A(ma,nb,mb,na) = M2_A(ma,nb,mb,na) & 
  +1.0d0 * t1_A(ma, na) * t1_A(nb, mb)

  ! ### Spin case: i_b, j_a, a_a, b_b ###

  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      do ja = i_ja, f_ja
        if (ja == ma) cycle 
        jb = ja + cc_nOa
        do ib = i_ib, f_ib
          if (ib == nb) cycle 
          ia = ib - cc_nOa
          M2_A(ib,ja,aa,bb) = M2_A(ib,ja,aa,bb) & 
          -1.0d0 * t1_A(ja, aa) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) & 
          -1.0d0 * t1_A(ja, aa) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) & 
          +1.0d0 * t1_A(ja, aa) * t1_A(ma, na) * t2_A(ia, nb, ba, mb) & 
          -1.0d0 * t1_A(ja, aa) * t1_A(nb, mb) * t2_A(ma, ia, ba, na) & 
          -1.0d0 * t1_A(ja, aa) * t1_A(ma, ba) * t1_A(ia, na) * t1_A(nb, mb) & 
          -1.0d0 * t1_A(ib, bb) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) & 
          -1.0d0 * t1_A(ib, bb) * t1_A(ma, na) * t2_A(nb, jb, ab, mb) & 
          -1.0d0 * t1_A(ib, bb) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) & 
          +1.0d0 * t1_A(ib, bb) * t1_A(nb, mb) * t2_A(ma, jb, na, ab) & 
          -1.0d0 * t1_A(ib, bb) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(jb, mb) & 
          +1.0d0 * t1_A(jb, ab) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) & 
          +1.0d0 * t1_A(jb, ab) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) & 
          -1.0d0 * t1_A(jb, ab) * t1_A(ma, na) * t2_A(ia, nb, ba, mb) & 
          +1.0d0 * t1_A(jb, ab) * t1_A(nb, mb) * t2_A(ma, ia, ba, na) & 
          +1.0d0 * t1_A(jb, ab) * t1_A(ma, ba) * t1_A(ia, na) * t1_A(nb, mb) & 
          +1.0d0 * t1_A(ia, ba) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) & 
          +1.0d0 * t1_A(ia, ba) * t1_A(ma, na) * t2_A(nb, jb, ab, mb) & 
          +1.0d0 * t1_A(ia, ba) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) & 
          -1.0d0 * t1_A(ia, ba) * t1_A(nb, mb) * t2_A(ma, jb, na, ab) & 
          +1.0d0 * t1_A(ia, ba) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(jb, mb)
        enddo
      enddo
    enddo
  enddo

  !! Deltas:((mb, bb))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    do ja = i_ja, f_ja
      if (ja == ma) cycle 
      jb = ja + cc_nOa
      do ib = i_ib, f_ib
        if (ib == nb) cycle 
        ia = ib - cc_nOa
        M2_A(ib,ja,aa,mb) = M2_A(ib,ja,aa,mb) & 
        -1.0d0 * t1_A(ja, aa) * t2_A(ia, nb, na, mb) & 
        -1.0d0 * t1_A(ja, aa) * t1_A(ia, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(ib, mb) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) & 
        -1.0d0 * t1_A(ib, mb) * t1_A(ma, na) * t2_A(nb, jb, ab, mb) & 
        -1.0d0 * t1_A(ib, mb) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) & 
        +1.0d0 * t1_A(ib, mb) * t1_A(nb, mb) * t2_A(ma, jb, na, ab) & 
        -1.0d0 * t1_A(ib, mb) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(jb, mb) & 
        +1.0d0 * t1_A(jb, ab) * t2_A(ia, nb, na, mb) & 
        +1.0d0 * t1_A(jb, ab) * t1_A(ia, na) * t1_A(nb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((nb, ib))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      do ja = i_ja, f_ja
        if (ja == ma) cycle 
        jb = ja + cc_nOa
        M2_A(nb,ja,aa,bb) = M2_A(nb,ja,aa,bb) & 
        +1.0d0 * t1_A(ja, aa) * t2_A(ma, nb, ba, mb) & 
        +1.0d0 * t1_A(ja, aa) * t1_A(ma, ba) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(nb, bb) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) & 
        -1.0d0 * t1_A(nb, bb) * t1_A(ma, na) * t2_A(nb, jb, ab, mb) & 
        -1.0d0 * t1_A(nb, bb) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) & 
        +1.0d0 * t1_A(nb, bb) * t1_A(nb, mb) * t2_A(ma, jb, na, ab) & 
        -1.0d0 * t1_A(nb, bb) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(jb, mb) & 
        -1.0d0 * t1_A(jb, ab) * t2_A(ma, nb, ba, mb) & 
        -1.0d0 * t1_A(jb, ab) * t1_A(ma, ba) * t1_A(nb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((na, aa))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do ja = i_ja, f_ja
      if (ja == ma) cycle 
      jb = ja + cc_nOa
      do ib = i_ib, f_ib
        if (ib == nb) cycle 
        ia = ib - cc_nOa
        M2_A(ib,ja,na,bb) = M2_A(ib,ja,na,bb) & 
        -1.0d0 * t1_A(ja, na) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) & 
        -1.0d0 * t1_A(ja, na) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) & 
        +1.0d0 * t1_A(ja, na) * t1_A(ma, na) * t2_A(ia, nb, ba, mb) & 
        -1.0d0 * t1_A(ja, na) * t1_A(nb, mb) * t2_A(ma, ia, ba, na) & 
        -1.0d0 * t1_A(ja, na) * t1_A(ma, ba) * t1_A(ia, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(ib, bb) * t2_A(ma, jb, na, mb) & 
        -1.0d0 * t1_A(ib, bb) * t1_A(ma, na) * t1_A(jb, mb) & 
        +1.0d0 * t1_A(ia, ba) * t2_A(ma, jb, na, mb) & 
        +1.0d0 * t1_A(ia, ba) * t1_A(ma, na) * t1_A(jb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((ma, ja))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      do ib = i_ib, f_ib
        if (ib == nb) cycle 
        ia = ib - cc_nOa
        M2_A(ib,ma,aa,bb) = M2_A(ib,ma,aa,bb) & 
        -1.0d0 * t1_A(ma, aa) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) & 
        -1.0d0 * t1_A(ma, aa) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) & 
        +1.0d0 * t1_A(ma, aa) * t1_A(ma, na) * t2_A(ia, nb, ba, mb) & 
        -1.0d0 * t1_A(ma, aa) * t1_A(nb, mb) * t2_A(ma, ia, ba, na) & 
        -1.0d0 * t1_A(ma, aa) * t1_A(ma, ba) * t1_A(ia, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(ib, bb) * t2_A(ma, nb, na, ab) & 
        +1.0d0 * t1_A(ib, bb) * t1_A(nb, ab) * t1_A(ma, na) & 
        -1.0d0 * t1_A(ia, ba) * t2_A(ma, nb, na, ab) & 
        -1.0d0 * t1_A(ia, ba) * t1_A(nb, ab) * t1_A(ma, na)
      enddo
    enddo
  enddo

  !! Deltas:((nb, ib), (mb, bb))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    do ja = i_ja, f_ja
      if (ja == ma) cycle 
      jb = ja + cc_nOa
      M2_A(nb,ja,aa,mb) = M2_A(nb,ja,aa,mb) & 
      +1.0d0 * t1_A(ja, aa) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(ma, na) * t2_A(nb, jb, ab, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(nb, mb) * t2_A(ma, jb, na, ab) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(nb, ab) * t1_A(ma, na) * t1_A(jb, mb) & 
      -1.0d0 * t1_A(jb, ab) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((mb, bb), (na, aa))
  do ja = i_ja, f_ja
    if (ja == ma) cycle 
    jb = ja + cc_nOa
    do ib = i_ib, f_ib
      if (ib == nb) cycle 
      ia = ib - cc_nOa
      M2_A(ib,ja,na,mb) = M2_A(ib,ja,na,mb) & 
      -1.0d0 * t1_A(ja, na) * t2_A(ia, nb, na, mb) & 
      -1.0d0 * t1_A(ja, na) * t1_A(ia, na) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((nb, ib), (na, aa))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do ja = i_ja, f_ja
      if (ja == ma) cycle 
      jb = ja + cc_nOa
      M2_A(nb,ja,na,bb) = M2_A(nb,ja,na,bb) & 
      +1.0d0 * t1_A(ja, na) * t2_A(ma, nb, ba, mb) & 
      +1.0d0 * t1_A(ja, na) * t1_A(ma, ba) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(nb, bb) * t2_A(ma, jb, na, mb) & 
      -1.0d0 * t1_A(nb, bb) * t1_A(ma, na) * t1_A(jb, mb)
    enddo
  enddo

  !! Deltas:((ma, ja), (mb, bb))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    do ib = i_ib, f_ib
      if (ib == nb) cycle 
      ia = ib - cc_nOa
      M2_A(ib,ma,aa,mb) = M2_A(ib,ma,aa,mb) & 
      -1.0d0 * t1_A(ma, aa) * t2_A(ia, nb, na, mb) & 
      -1.0d0 * t1_A(ma, aa) * t1_A(ia, na) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(ib, mb) * t2_A(ma, nb, na, ab) & 
      +1.0d0 * t1_A(ib, mb) * t1_A(nb, ab) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((ma, ja), (nb, ib))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      M2_A(nb,ma,aa,bb) = M2_A(nb,ma,aa,bb) & 
      +1.0d0 * t1_A(ma, aa) * t2_A(ma, nb, ba, mb) & 
      +1.0d0 * t1_A(ma, aa) * t1_A(ma, ba) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((ma, ja), (na, aa))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do ib = i_ib, f_ib
      if (ib == nb) cycle 
      ia = ib - cc_nOa
      M2_A(ib,ma,na,bb) = M2_A(ib,ma,na,bb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(ma, na) * t2_A(ia, nb, ba, mb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(nb, mb) * t2_A(ma, ia, ba, na) & 
      -1.0d0 * t1_A(ma, na) * t1_A(ma, ba) * t1_A(ia, na) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(ib, bb) * t1_A(ma, na) & 
      -1.0d0 * t1_A(ia, ba) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((na, aa), (mb, bb))
  do ja = i_ja, f_ja
    if (ja == ma) cycle 
    jb = ja + cc_nOa
    do ib = i_ib, f_ib
      if (ib == nb) cycle 
      ia = ib - cc_nOa
      M2_A(ib,ja,na,mb) = M2_A(ib,ja,na,mb) & 
      -1.0d0 * t1_A(ib, mb) * t2_A(ma, jb, na, mb) & 
      -1.0d0 * t1_A(ib, mb) * t1_A(ma, na) * t1_A(jb, mb)
    enddo
  enddo

  !! Deltas:((nb, ib), (ma, ja))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do aa = i_aa, f_aa
      if (aa == na) cycle 
      ab = aa + cc_nVa
      M2_A(nb,ma,aa,bb) = M2_A(nb,ma,aa,bb) & 
      +1.0d0 * t1_A(nb, bb) * t2_A(ma, nb, na, ab) & 
      +1.0d0 * t1_A(nb, bb) * t1_A(nb, ab) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((nb, ib), (na, aa), (mb, bb))
  do ja = i_ja, f_ja
    if (ja == ma) cycle 
    jb = ja + cc_nOa
    M2_A(nb,ja,na,mb) = M2_A(nb,ja,na,mb) & 
    +1.0d0 * t1_A(ja, na) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ja), (mb, bb), (nb, ib))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    M2_A(nb,ma,aa,mb) = M2_A(nb,ma,aa,mb) & 
    +1.0d0 * t1_A(ma, aa) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ja), (na, aa), (mb, bb))
  do ib = i_ib, f_ib
    if (ib == nb) cycle 
    ia = ib - cc_nOa
    M2_A(ib,ma,na,mb) = M2_A(ib,ma,na,mb) & 
    -1.0d0 * t1_A(ma, na) * t2_A(ia, nb, na, mb) & 
    -1.0d0 * t1_A(ma, na) * t1_A(ia, na) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ja), (na, aa), (nb, ib))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    M2_A(nb,ma,na,bb) = M2_A(nb,ma,na,bb) & 
    +1.0d0 * t1_A(ma, na) * t2_A(ma, nb, ba, mb) & 
    +1.0d0 * t1_A(ma, na) * t1_A(ma, ba) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ja), (mb, bb), (na, aa))
  do ib = i_ib, f_ib
    if (ib == nb) cycle 
    ia = ib - cc_nOa
    M2_A(ib,ma,na,mb) = M2_A(ib,ma,na,mb) & 
    +1.0d0 * t1_A(ib, mb) * t1_A(ma, na)
  enddo

  !! Deltas:((nb, ib), (na, aa), (ma, ja))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    M2_A(nb,ma,na,bb) = M2_A(nb,ma,na,bb) & 
    +1.0d0 * t1_A(nb, bb) * t1_A(ma, na)
  enddo

  !! Deltas:((nb, ib), (mb, bb), (na, aa))
  do ja = i_ja, f_ja
    if (ja == ma) cycle 
    jb = ja + cc_nOa
    M2_A(nb,ja,na,mb) = M2_A(nb,ja,na,mb) & 
    -1.0d0 * t1_A(nb, mb) * t2_A(ma, jb, na, mb) & 
    -1.0d0 * t1_A(nb, mb) * t1_A(ma, na) * t1_A(jb, mb)
  enddo

  !! Deltas:((nb, ib), (mb, bb), (ma, ja))
  do aa = i_aa, f_aa
    if (aa == na) cycle 
    ab = aa + cc_nVa
    M2_A(nb,ma,aa,mb) = M2_A(nb,ma,aa,mb) & 
    +1.0d0 * t1_A(nb, mb) * t2_A(ma, nb, na, ab) & 
    +1.0d0 * t1_A(nb, mb) * t1_A(nb, ab) * t1_A(ma, na)
  enddo

  !! Deltas:((ma, ja), (na, aa), (mb, bb), (nb, ib))
  M2_A(nb,ma,na,mb) = M2_A(nb,ma,na,mb) & 
  +1.0d0 * t1_A(ma, na) * t1_A(nb, mb)

  !! Deltas:((nb, ib), (mb, bb), (na, aa), (ma, ja))
  M2_A(nb,ma,na,mb) = M2_A(nb,ma,na,mb) & 
  +1.0d0 * t1_A(nb, mb) * t1_A(ma, na)

  ! ### Spin case: i_b, j_a, a_b, b_a ###

  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      do ja = i_ja, f_ja
        if (ja == ma) cycle 
        jb = ja + cc_nOa
        do ib = i_ib, f_ib
          if (ib == nb) cycle 
          ia = ib - cc_nOa
          M2_A(ib,ja,ab,ba) = M2_A(ib,ja,ab,ba) & 
          +1.0d0 * t1_A(ib, ab) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) & 
          +1.0d0 * t1_A(ib, ab) * t1_A(ma, na) * t2_A(nb, jb, bb, mb) & 
          +1.0d0 * t1_A(ib, ab) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) & 
          -1.0d0 * t1_A(ib, ab) * t1_A(nb, mb) * t2_A(ma, jb, na, bb) & 
          +1.0d0 * t1_A(ib, ab) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(jb, mb) & 
          +1.0d0 * t1_A(ja, ba) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) & 
          +1.0d0 * t1_A(ja, ba) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) & 
          -1.0d0 * t1_A(ja, ba) * t1_A(ma, na) * t2_A(ia, nb, aa, mb) & 
          +1.0d0 * t1_A(ja, ba) * t1_A(nb, mb) * t2_A(ma, ia, aa, na) & 
          +1.0d0 * t1_A(ja, ba) * t1_A(ma, aa) * t1_A(ia, na) * t1_A(nb, mb) & 
          -1.0d0 * t1_A(ia, aa) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) & 
          -1.0d0 * t1_A(ia, aa) * t1_A(ma, na) * t2_A(nb, jb, bb, mb) & 
          -1.0d0 * t1_A(ia, aa) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) & 
          +1.0d0 * t1_A(ia, aa) * t1_A(nb, mb) * t2_A(ma, jb, na, bb) & 
          -1.0d0 * t1_A(ia, aa) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(jb, mb) & 
          -1.0d0 * t1_A(jb, bb) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) & 
          -1.0d0 * t1_A(jb, bb) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) & 
          +1.0d0 * t1_A(jb, bb) * t1_A(ma, na) * t2_A(ia, nb, aa, mb) & 
          -1.0d0 * t1_A(jb, bb) * t1_A(nb, mb) * t2_A(ma, ia, aa, na) & 
          -1.0d0 * t1_A(jb, bb) * t1_A(ma, aa) * t1_A(ia, na) * t1_A(nb, mb)
        enddo
      enddo
    enddo
  enddo

  !! Deltas:((na, ba))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    do ja = i_ja, f_ja
      if (ja == ma) cycle 
      jb = ja + cc_nOa
      do ib = i_ib, f_ib
        if (ib == nb) cycle 
        ia = ib - cc_nOa
        M2_A(ib,ja,ab,na) = M2_A(ib,ja,ab,na) & 
        +1.0d0 * t1_A(ib, ab) * t2_A(ma, jb, na, mb) & 
        +1.0d0 * t1_A(ib, ab) * t1_A(ma, na) * t1_A(jb, mb) & 
        +1.0d0 * t1_A(ja, na) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) & 
        +1.0d0 * t1_A(ja, na) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) & 
        -1.0d0 * t1_A(ja, na) * t1_A(ma, na) * t2_A(ia, nb, aa, mb) & 
        +1.0d0 * t1_A(ja, na) * t1_A(nb, mb) * t2_A(ma, ia, aa, na) & 
        +1.0d0 * t1_A(ja, na) * t1_A(ma, aa) * t1_A(ia, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(ia, aa) * t2_A(ma, jb, na, mb) & 
        -1.0d0 * t1_A(ia, aa) * t1_A(ma, na) * t1_A(jb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((ma, ja))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      do ib = i_ib, f_ib
        if (ib == nb) cycle 
        ia = ib - cc_nOa
        M2_A(ib,ma,ab,ba) = M2_A(ib,ma,ab,ba) & 
        -1.0d0 * t1_A(ib, ab) * t2_A(ma, nb, na, bb) & 
        -1.0d0 * t1_A(ib, ab) * t1_A(nb, bb) * t1_A(ma, na) & 
        +1.0d0 * t1_A(ma, ba) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) & 
        +1.0d0 * t1_A(ma, ba) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) & 
        -1.0d0 * t1_A(ma, ba) * t1_A(ma, na) * t2_A(ia, nb, aa, mb) & 
        +1.0d0 * t1_A(ma, ba) * t1_A(nb, mb) * t2_A(ma, ia, aa, na) & 
        +1.0d0 * t1_A(ma, ba) * t1_A(ma, aa) * t1_A(ia, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(ia, aa) * t2_A(ma, nb, na, bb) & 
        +1.0d0 * t1_A(ia, aa) * t1_A(nb, bb) * t1_A(ma, na)
      enddo
    enddo
  enddo

  !! Deltas:((mb, ab))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ja = i_ja, f_ja
      if (ja == ma) cycle 
      jb = ja + cc_nOa
      do ib = i_ib, f_ib
        if (ib == nb) cycle 
        ia = ib - cc_nOa
        M2_A(ib,ja,mb,ba) = M2_A(ib,ja,mb,ba) & 
        +1.0d0 * t1_A(ib, mb) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) & 
        +1.0d0 * t1_A(ib, mb) * t1_A(ma, na) * t2_A(nb, jb, bb, mb) & 
        +1.0d0 * t1_A(ib, mb) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) & 
        -1.0d0 * t1_A(ib, mb) * t1_A(nb, mb) * t2_A(ma, jb, na, bb) & 
        +1.0d0 * t1_A(ib, mb) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(jb, mb) & 
        +1.0d0 * t1_A(ja, ba) * t2_A(ia, nb, na, mb) & 
        +1.0d0 * t1_A(ja, ba) * t1_A(ia, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(jb, bb) * t2_A(ia, nb, na, mb) & 
        -1.0d0 * t1_A(jb, bb) * t1_A(ia, na) * t1_A(nb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((nb, ib))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      do ja = i_ja, f_ja
        if (ja == ma) cycle 
        jb = ja + cc_nOa
        M2_A(nb,ja,ab,ba) = M2_A(nb,ja,ab,ba) & 
        +1.0d0 * t1_A(nb, ab) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) & 
        +1.0d0 * t1_A(nb, ab) * t1_A(ma, na) * t2_A(nb, jb, bb, mb) & 
        +1.0d0 * t1_A(nb, ab) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) & 
        -1.0d0 * t1_A(nb, ab) * t1_A(nb, mb) * t2_A(ma, jb, na, bb) & 
        +1.0d0 * t1_A(nb, ab) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(jb, mb) & 
        -1.0d0 * t1_A(ja, ba) * t2_A(ma, nb, aa, mb) & 
        -1.0d0 * t1_A(ja, ba) * t1_A(ma, aa) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(jb, bb) * t2_A(ma, nb, aa, mb) & 
        +1.0d0 * t1_A(jb, bb) * t1_A(ma, aa) * t1_A(nb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((ma, ja), (na, ba))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    do ib = i_ib, f_ib
      if (ib == nb) cycle 
      ia = ib - cc_nOa
      M2_A(ib,ma,ab,na) = M2_A(ib,ma,ab,na) & 
      -1.0d0 * t1_A(ib, ab) * t1_A(ma, na) & 
      +1.0d0 * t1_A(ma, na) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) & 
      -1.0d0 * t1_A(ma, na) * t1_A(ma, na) * t2_A(ia, nb, aa, mb) & 
      +1.0d0 * t1_A(ma, na) * t1_A(nb, mb) * t2_A(ma, ia, aa, na) & 
      +1.0d0 * t1_A(ma, na) * t1_A(ma, aa) * t1_A(ia, na) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(ia, aa) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((na, ba), (mb, ab))
  do ja = i_ja, f_ja
    if (ja == ma) cycle 
    jb = ja + cc_nOa
    do ib = i_ib, f_ib
      if (ib == nb) cycle 
      ia = ib - cc_nOa
      M2_A(ib,ja,mb,na) = M2_A(ib,ja,mb,na) & 
      +1.0d0 * t1_A(ib, mb) * t2_A(ma, jb, na, mb) & 
      +1.0d0 * t1_A(ib, mb) * t1_A(ma, na) * t1_A(jb, mb)
    enddo
  enddo

  !! Deltas:((ma, ja), (mb, ab))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ib = i_ib, f_ib
      if (ib == nb) cycle 
      ia = ib - cc_nOa
      M2_A(ib,ma,mb,ba) = M2_A(ib,ma,mb,ba) & 
      -1.0d0 * t1_A(ib, mb) * t2_A(ma, nb, na, bb) & 
      -1.0d0 * t1_A(ib, mb) * t1_A(nb, bb) * t1_A(ma, na) & 
      +1.0d0 * t1_A(ma, ba) * t2_A(ia, nb, na, mb) & 
      +1.0d0 * t1_A(ma, ba) * t1_A(ia, na) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((nb, ib), (na, ba))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    do ja = i_ja, f_ja
      if (ja == ma) cycle 
      jb = ja + cc_nOa
      M2_A(nb,ja,ab,na) = M2_A(nb,ja,ab,na) & 
      +1.0d0 * t1_A(nb, ab) * t2_A(ma, jb, na, mb) & 
      +1.0d0 * t1_A(nb, ab) * t1_A(ma, na) * t1_A(jb, mb) & 
      -1.0d0 * t1_A(ja, na) * t2_A(ma, nb, aa, mb) & 
      -1.0d0 * t1_A(ja, na) * t1_A(ma, aa) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((nb, ib), (ma, ja))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      M2_A(nb,ma,ab,ba) = M2_A(nb,ma,ab,ba) & 
      -1.0d0 * t1_A(nb, ab) * t2_A(ma, nb, na, bb) & 
      -1.0d0 * t1_A(nb, ab) * t1_A(nb, bb) * t1_A(ma, na)
    enddo
  enddo

  !! Deltas:((nb, ib), (mb, ab))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ja = i_ja, f_ja
      if (ja == ma) cycle 
      jb = ja + cc_nOa
      M2_A(nb,ja,mb,ba) = M2_A(nb,ja,mb,ba) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(ma, na) * t2_A(nb, jb, bb, mb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(nb, mb) * t2_A(ma, jb, na, bb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(nb, bb) * t1_A(ma, na) * t1_A(jb, mb) & 
      -1.0d0 * t1_A(ja, ba) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(jb, bb) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((mb, ab), (na, ba))
  do ja = i_ja, f_ja
    if (ja == ma) cycle 
    jb = ja + cc_nOa
    do ib = i_ib, f_ib
      if (ib == nb) cycle 
      ia = ib - cc_nOa
      M2_A(ib,ja,mb,na) = M2_A(ib,ja,mb,na) & 
      +1.0d0 * t1_A(ja, na) * t2_A(ia, nb, na, mb) & 
      +1.0d0 * t1_A(ja, na) * t1_A(ia, na) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((ma, ja), (nb, ib))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      M2_A(nb,ma,ab,ba) = M2_A(nb,ma,ab,ba) & 
      -1.0d0 * t1_A(ma, ba) * t2_A(ma, nb, aa, mb) & 
      -1.0d0 * t1_A(ma, ba) * t1_A(ma, aa) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((ma, ja), (mb, ab), (na, ba))
  do ib = i_ib, f_ib
    if (ib == nb) cycle 
    ia = ib - cc_nOa
    M2_A(ib,ma,mb,na) = M2_A(ib,ma,mb,na) & 
    -1.0d0 * t1_A(ib, mb) * t1_A(ma, na)
  enddo

  !! Deltas:((nb, ib), (na, ba), (ma, ja))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    M2_A(nb,ma,ab,na) = M2_A(nb,ma,ab,na) & 
    -1.0d0 * t1_A(nb, ab) * t1_A(ma, na)
  enddo

  !! Deltas:((nb, ib), (mb, ab), (na, ba))
  do ja = i_ja, f_ja
    if (ja == ma) cycle 
    jb = ja + cc_nOa
    M2_A(nb,ja,mb,na) = M2_A(nb,ja,mb,na) & 
    +1.0d0 * t1_A(nb, mb) * t2_A(ma, jb, na, mb) & 
    +1.0d0 * t1_A(nb, mb) * t1_A(ma, na) * t1_A(jb, mb)
  enddo

  !! Deltas:((nb, ib), (mb, ab), (ma, ja))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    M2_A(nb,ma,mb,ba) = M2_A(nb,ma,mb,ba) & 
    -1.0d0 * t1_A(nb, mb) * t2_A(ma, nb, na, bb) & 
    -1.0d0 * t1_A(nb, mb) * t1_A(nb, bb) * t1_A(ma, na)
  enddo

  !! Deltas:((nb, ib), (na, ba), (mb, ab))
  do ja = i_ja, f_ja
    if (ja == ma) cycle 
    jb = ja + cc_nOa
    M2_A(nb,ja,mb,na) = M2_A(nb,ja,mb,na) & 
    -1.0d0 * t1_A(ja, na) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ja), (mb, ab), (nb, ib))
  do ba = i_ba, f_ba
    if (ba == na) cycle 
    bb = ba + cc_nVa
    M2_A(nb,ma,mb,ba) = M2_A(nb,ma,mb,ba) & 
    -1.0d0 * t1_A(ma, ba) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ja), (na, ba), (mb, ab))
  do ib = i_ib, f_ib
    if (ib == nb) cycle 
    ia = ib - cc_nOa
    M2_A(ib,ma,mb,na) = M2_A(ib,ma,mb,na) & 
    +1.0d0 * t1_A(ma, na) * t2_A(ia, nb, na, mb) & 
    +1.0d0 * t1_A(ma, na) * t1_A(ia, na) * t1_A(nb, mb)
  enddo

  !! Deltas:((ma, ja), (na, ba), (nb, ib))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    M2_A(nb,ma,ab,na) = M2_A(nb,ma,ab,na) & 
    -1.0d0 * t1_A(ma, na) * t2_A(ma, nb, aa, mb) & 
    -1.0d0 * t1_A(ma, na) * t1_A(ma, aa) * t1_A(nb, mb)
  enddo

  !! Deltas:((nb, ib), (mb, ab), (na, ba), (ma, ja))
  M2_A(nb,ma,mb,na) = M2_A(nb,ma,mb,na) & 
  -1.0d0 * t1_A(nb, mb) * t1_A(ma, na)

  !! Deltas:((ma, ja), (na, ba), (mb, ab), (nb, ib))
  M2_A(nb,ma,mb,na) = M2_A(nb,ma,mb,na) & 
  -1.0d0 * t1_A(ma, na) * t1_A(nb, mb)

  ! ### Spin case: i_b, j_b, a_b, b_b ###

  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      do jb = i_jb, f_jb
        if (jb == nb) cycle 
        ja = jb - cc_nOa
        do ib = i_ib, f_ib
          if (ib == nb) cycle 
          ia = ib - cc_nOa
          M2_A(ib,jb,ab,bb) = M2_A(ib,jb,ab,bb) & 
          +1.0d0 * t1_A(ib, ab) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) & 
          +1.0d0 * t1_A(ib, ab) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) & 
          -1.0d0 * t1_A(ib, ab) * t1_A(ma, na) * t2_A(ja, nb, ba, mb) & 
          +1.0d0 * t1_A(ib, ab) * t1_A(nb, mb) * t2_A(ma, ja, ba, na) & 
          +1.0d0 * t1_A(ib, ab) * t1_A(ma, ba) * t1_A(ja, na) * t1_A(nb, mb) & 
          -1.0d0 * t1_A(jb, ab) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) & 
          -1.0d0 * t1_A(jb, ab) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) & 
          +1.0d0 * t1_A(jb, ab) * t1_A(ma, na) * t2_A(ia, nb, ba, mb) & 
          -1.0d0 * t1_A(jb, ab) * t1_A(nb, mb) * t2_A(ma, ia, ba, na) & 
          -1.0d0 * t1_A(jb, ab) * t1_A(ma, ba) * t1_A(ia, na) * t1_A(nb, mb) & 
          -1.0d0 * t1_A(ib, bb) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) & 
          -1.0d0 * t1_A(ib, bb) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) & 
          +1.0d0 * t1_A(ib, bb) * t1_A(ma, na) * t2_A(ja, nb, aa, mb) & 
          -1.0d0 * t1_A(ib, bb) * t1_A(nb, mb) * t2_A(ma, ja, aa, na) & 
          -1.0d0 * t1_A(ib, bb) * t1_A(ma, aa) * t1_A(ja, na) * t1_A(nb, mb) & 
          +1.0d0 * t1_A(jb, bb) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) & 
          +1.0d0 * t1_A(jb, bb) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) & 
          -1.0d0 * t1_A(jb, bb) * t1_A(ma, na) * t2_A(ia, nb, aa, mb) & 
          +1.0d0 * t1_A(jb, bb) * t1_A(nb, mb) * t2_A(ma, ia, aa, na) & 
          +1.0d0 * t1_A(jb, bb) * t1_A(ma, aa) * t1_A(ia, na) * t1_A(nb, mb) & 
          -1.0d0 * t1_A(ia, aa) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) & 
          -1.0d0 * t1_A(ia, aa) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) & 
          +1.0d0 * t1_A(ia, aa) * t1_A(ma, na) * t2_A(ja, nb, ba, mb) & 
          -1.0d0 * t1_A(ia, aa) * t1_A(nb, mb) * t2_A(ma, ja, ba, na) & 
          -1.0d0 * t1_A(ia, aa) * t1_A(ma, ba) * t1_A(ja, na) * t1_A(nb, mb) & 
          +1.0d0 * t1_A(ja, aa) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) & 
          +1.0d0 * t1_A(ja, aa) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) & 
          -1.0d0 * t1_A(ja, aa) * t1_A(ma, na) * t2_A(ia, nb, ba, mb) & 
          +1.0d0 * t1_A(ja, aa) * t1_A(nb, mb) * t2_A(ma, ia, ba, na) & 
          +1.0d0 * t1_A(ja, aa) * t1_A(ma, ba) * t1_A(ia, na) * t1_A(nb, mb) & 
          +1.0d0 * t1_A(ia, ba) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) & 
          +1.0d0 * t1_A(ia, ba) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) & 
          -1.0d0 * t1_A(ia, ba) * t1_A(ma, na) * t2_A(ja, nb, aa, mb) & 
          +1.0d0 * t1_A(ia, ba) * t1_A(nb, mb) * t2_A(ma, ja, aa, na) & 
          +1.0d0 * t1_A(ia, ba) * t1_A(ma, aa) * t1_A(ja, na) * t1_A(nb, mb) & 
          -1.0d0 * t1_A(ja, ba) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) & 
          -1.0d0 * t1_A(ja, ba) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) & 
          +1.0d0 * t1_A(ja, ba) * t1_A(ma, na) * t2_A(ia, nb, aa, mb) & 
          -1.0d0 * t1_A(ja, ba) * t1_A(nb, mb) * t2_A(ma, ia, aa, na) & 
          -1.0d0 * t1_A(ja, ba) * t1_A(ma, aa) * t1_A(ia, na) * t1_A(nb, mb)
        enddo
      enddo
    enddo
  enddo

  !! Deltas:((mb, bb))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    do jb = i_jb, f_jb
      if (jb == nb) cycle 
      ja = jb - cc_nOa
      do ib = i_ib, f_ib
        if (ib == nb) cycle 
        ia = ib - cc_nOa
        M2_A(ib,jb,ab,mb) = M2_A(ib,jb,ab,mb) & 
        +1.0d0 * t1_A(ib, ab) * t2_A(ja, nb, na, mb) & 
        +1.0d0 * t1_A(ib, ab) * t1_A(ja, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(jb, ab) * t2_A(ia, nb, na, mb) & 
        -1.0d0 * t1_A(jb, ab) * t1_A(ia, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(ib, mb) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) & 
        -1.0d0 * t1_A(ib, mb) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) & 
        +1.0d0 * t1_A(ib, mb) * t1_A(ma, na) * t2_A(ja, nb, aa, mb) & 
        -1.0d0 * t1_A(ib, mb) * t1_A(nb, mb) * t2_A(ma, ja, aa, na) & 
        -1.0d0 * t1_A(ib, mb) * t1_A(ma, aa) * t1_A(ja, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(jb, mb) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) & 
        +1.0d0 * t1_A(jb, mb) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) & 
        -1.0d0 * t1_A(jb, mb) * t1_A(ma, na) * t2_A(ia, nb, aa, mb) & 
        +1.0d0 * t1_A(jb, mb) * t1_A(nb, mb) * t2_A(ma, ia, aa, na) & 
        +1.0d0 * t1_A(jb, mb) * t1_A(ma, aa) * t1_A(ia, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(ia, aa) * t2_A(ja, nb, na, mb) & 
        -1.0d0 * t1_A(ia, aa) * t1_A(ja, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(ja, aa) * t2_A(ia, nb, na, mb) & 
        +1.0d0 * t1_A(ja, aa) * t1_A(ia, na) * t1_A(nb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((nb, jb))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      do ib = i_ib, f_ib
        if (ib == nb) cycle 
        ia = ib - cc_nOa
        M2_A(ib,nb,ab,bb) = M2_A(ib,nb,ab,bb) & 
        -1.0d0 * t1_A(ib, ab) * t2_A(ma, nb, ba, mb) & 
        -1.0d0 * t1_A(ib, ab) * t1_A(ma, ba) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(nb, ab) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) & 
        -1.0d0 * t1_A(nb, ab) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) & 
        +1.0d0 * t1_A(nb, ab) * t1_A(ma, na) * t2_A(ia, nb, ba, mb) & 
        -1.0d0 * t1_A(nb, ab) * t1_A(nb, mb) * t2_A(ma, ia, ba, na) & 
        -1.0d0 * t1_A(nb, ab) * t1_A(ma, ba) * t1_A(ia, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(ib, bb) * t2_A(ma, nb, aa, mb) & 
        +1.0d0 * t1_A(ib, bb) * t1_A(ma, aa) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(nb, bb) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) & 
        +1.0d0 * t1_A(nb, bb) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) & 
        -1.0d0 * t1_A(nb, bb) * t1_A(ma, na) * t2_A(ia, nb, aa, mb) & 
        +1.0d0 * t1_A(nb, bb) * t1_A(nb, mb) * t2_A(ma, ia, aa, na) & 
        +1.0d0 * t1_A(nb, bb) * t1_A(ma, aa) * t1_A(ia, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(ia, aa) * t2_A(ma, nb, ba, mb) & 
        +1.0d0 * t1_A(ia, aa) * t1_A(ma, ba) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(ia, ba) * t2_A(ma, nb, aa, mb) & 
        -1.0d0 * t1_A(ia, ba) * t1_A(ma, aa) * t1_A(nb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((mb, ab))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do jb = i_jb, f_jb
      if (jb == nb) cycle 
      ja = jb - cc_nOa
      do ib = i_ib, f_ib
        if (ib == nb) cycle 
        ia = ib - cc_nOa
        M2_A(ib,jb,mb,bb) = M2_A(ib,jb,mb,bb) & 
        +1.0d0 * t1_A(ib, mb) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) & 
        +1.0d0 * t1_A(ib, mb) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) & 
        -1.0d0 * t1_A(ib, mb) * t1_A(ma, na) * t2_A(ja, nb, ba, mb) & 
        +1.0d0 * t1_A(ib, mb) * t1_A(nb, mb) * t2_A(ma, ja, ba, na) & 
        +1.0d0 * t1_A(ib, mb) * t1_A(ma, ba) * t1_A(ja, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(jb, mb) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) & 
        -1.0d0 * t1_A(jb, mb) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) & 
        +1.0d0 * t1_A(jb, mb) * t1_A(ma, na) * t2_A(ia, nb, ba, mb) & 
        -1.0d0 * t1_A(jb, mb) * t1_A(nb, mb) * t2_A(ma, ia, ba, na) & 
        -1.0d0 * t1_A(jb, mb) * t1_A(ma, ba) * t1_A(ia, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(ib, bb) * t2_A(ja, nb, na, mb) & 
        -1.0d0 * t1_A(ib, bb) * t1_A(ja, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(jb, bb) * t2_A(ia, nb, na, mb) & 
        +1.0d0 * t1_A(jb, bb) * t1_A(ia, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(ia, ba) * t2_A(ja, nb, na, mb) & 
        +1.0d0 * t1_A(ia, ba) * t1_A(ja, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(ja, ba) * t2_A(ia, nb, na, mb) & 
        -1.0d0 * t1_A(ja, ba) * t1_A(ia, na) * t1_A(nb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((nb, ib))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do ab = i_ab, f_ab
      if (ab == mb) cycle 
      aa = ab - cc_nVa
      do jb = i_jb, f_jb
        if (jb == nb) cycle 
        ja = jb - cc_nOa
        M2_A(nb,jb,ab,bb) = M2_A(nb,jb,ab,bb) & 
        +1.0d0 * t1_A(nb, ab) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) & 
        +1.0d0 * t1_A(nb, ab) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) & 
        -1.0d0 * t1_A(nb, ab) * t1_A(ma, na) * t2_A(ja, nb, ba, mb) & 
        +1.0d0 * t1_A(nb, ab) * t1_A(nb, mb) * t2_A(ma, ja, ba, na) & 
        +1.0d0 * t1_A(nb, ab) * t1_A(ma, ba) * t1_A(ja, na) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(jb, ab) * t2_A(ma, nb, ba, mb) & 
        +1.0d0 * t1_A(jb, ab) * t1_A(ma, ba) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(nb, bb) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) & 
        -1.0d0 * t1_A(nb, bb) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) & 
        +1.0d0 * t1_A(nb, bb) * t1_A(ma, na) * t2_A(ja, nb, aa, mb) & 
        -1.0d0 * t1_A(nb, bb) * t1_A(nb, mb) * t2_A(ma, ja, aa, na) & 
        -1.0d0 * t1_A(nb, bb) * t1_A(ma, aa) * t1_A(ja, na) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(jb, bb) * t2_A(ma, nb, aa, mb) & 
        -1.0d0 * t1_A(jb, bb) * t1_A(ma, aa) * t1_A(nb, mb) & 
        -1.0d0 * t1_A(ja, aa) * t2_A(ma, nb, ba, mb) & 
        -1.0d0 * t1_A(ja, aa) * t1_A(ma, ba) * t1_A(nb, mb) & 
        +1.0d0 * t1_A(ja, ba) * t2_A(ma, nb, aa, mb) & 
        +1.0d0 * t1_A(ja, ba) * t1_A(ma, aa) * t1_A(nb, mb)
      enddo
    enddo
  enddo

  !! Deltas:((nb, jb), (mb, bb))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    do ib = i_ib, f_ib
      if (ib == nb) cycle 
      ia = ib - cc_nOa
      M2_A(ib,nb,ab,mb) = M2_A(ib,nb,ab,mb) & 
      -1.0d0 * t1_A(ib, ab) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(nb, ab) * t2_A(ia, nb, na, mb) & 
      -1.0d0 * t1_A(nb, ab) * t1_A(ia, na) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(ib, mb) * t2_A(ma, nb, aa, mb) & 
      +1.0d0 * t1_A(ib, mb) * t1_A(ma, aa) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(ma, na) * t2_A(ia, nb, aa, mb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(nb, mb) * t2_A(ma, ia, aa, na) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(ma, aa) * t1_A(ia, na) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(ia, aa) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((nb, jb), (mb, ab))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do ib = i_ib, f_ib
      if (ib == nb) cycle 
      ia = ib - cc_nOa
      M2_A(ib,nb,mb,bb) = M2_A(ib,nb,mb,bb) & 
      -1.0d0 * t1_A(ib, mb) * t2_A(ma, nb, ba, mb) & 
      -1.0d0 * t1_A(ib, mb) * t1_A(ma, ba) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(ma, na) * t2_A(ia, nb, ba, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(nb, mb) * t2_A(ma, ia, ba, na) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(ma, ba) * t1_A(ia, na) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(ib, bb) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(nb, bb) * t2_A(ia, nb, na, mb) & 
      +1.0d0 * t1_A(nb, bb) * t1_A(ia, na) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(ia, ba) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((nb, ib), (mb, bb))
  do ab = i_ab, f_ab
    if (ab == mb) cycle 
    aa = ab - cc_nVa
    do jb = i_jb, f_jb
      if (jb == nb) cycle 
      ja = jb - cc_nOa
      M2_A(nb,jb,ab,mb) = M2_A(nb,jb,ab,mb) & 
      +1.0d0 * t1_A(nb, ab) * t2_A(ja, nb, na, mb) & 
      +1.0d0 * t1_A(nb, ab) * t1_A(ja, na) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(jb, ab) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(ma, na) * t2_A(ja, nb, aa, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(nb, mb) * t2_A(ma, ja, aa, na) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(ma, aa) * t1_A(ja, na) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(jb, mb) * t2_A(ma, nb, aa, mb) & 
      -1.0d0 * t1_A(jb, mb) * t1_A(ma, aa) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(ja, aa) * t1_A(nb, mb)
    enddo
  enddo

  !! Deltas:((nb, ib), (mb, ab))
  do bb = i_bb, f_bb
    if (bb == mb) cycle 
    ba = bb - cc_nVa
    do jb = i_jb, f_jb
      if (jb == nb) cycle 
      ja = jb - cc_nOa
      M2_A(nb,jb,mb,bb) = M2_A(nb,jb,mb,bb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) & 
      -1.0d0 * t1_A(nb, mb) * t1_A(ma, na) * t2_A(ja, nb, ba, mb) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(nb, mb) * t2_A(ma, ja, ba, na) & 
      +1.0d0 * t1_A(nb, mb) * t1_A(ma, ba) * t1_A(ja, na) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(jb, mb) * t2_A(ma, nb, ba, mb) & 
      +1.0d0 * t1_A(jb, mb) * t1_A(ma, ba) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(nb, bb) * t2_A(ja, nb, na, mb) & 
      -1.0d0 * t1_A(nb, bb) * t1_A(ja, na) * t1_A(nb, mb) & 
      -1.0d0 * t1_A(jb, bb) * t1_A(nb, mb) & 
      +1.0d0 * t1_A(ja, ba) * t1_A(nb, mb)
    enddo
  enddo


end
