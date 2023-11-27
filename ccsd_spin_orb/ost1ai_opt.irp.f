subroutine ost1ai_opt(nO,nV,det,t1_A,t2_A,M1_A,M2_A)

  implicit none

  integer, intent(in)           :: nO,nV
  integer(bit_kind), intent(in) :: det(N_int,2)
  double precision, intent(in)  :: t1_A(nO,nV), t2_A(nO,nO,nV,nV)
  
  double precision, intent(out) :: M1_A(nO,nV),M2_A(nO,nO,nV,nV)

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

  !print*,'ia',i_ia,f_ia
  !print*,'ib',i_ib,f_ib
  !print*,'aa',i_aa,f_aa
  !print*,'ab',i_ab,f_ab
  
  ! ### Spin case: i_a, a_a ###

  !! Deltas:((na, aa))
  do ia = i_ia, f_ia
    if (ia == ma) cycle
    ib = ia + cc_nOa
    M1_A(ia,na) = M1_A(ia,na) &
    -1.0d0 * t2_A(ma, ib, na, mb)
  enddo

  ! ### Spin case: i_b, a_b ###

  !! Deltas:((mb, ab))
  do ib = i_ib, f_ib
    if (ib == nb) cycle
    ia = ib - cc_nOa
    M1_A(ib,mb) = M1_A(ib,mb) &
    -1.0d0 * t2_A(ia, nb, na, mb)
  enddo

  ! ### Spin case: i_a, j_a, a_a, b_a ###

  !! Deltas:((na, aa))
      do ba = i_ba, f_ba
!        if (ba == na) cycle
        bb = ba + cc_nVa
    do ja = i_ja, f_ja
!      if (ja == ma) cycle
      jb = ja + cc_nOa
  do ia = i_ia, f_ia
!    if (ia == ma) cycle
    ib = ia + cc_nOa
        M2_A(ia,ja,na,ba) = M2_A(ia,ja,na,ba) &
        !+1.0d0 * t1_A(ia, ba) * t2_A(ma, jb, na, mb) &
        !-1.0d0 * t1_A(ja, ba) * t2_A(ma, ib, na, mb) &
        !-1.0d0 * t1_A(ib, bb) * t2_A(ma, jb, na, mb) &
        !+1.0d0 * t1_A(jb, bb) * t2_A(ma, ib, na, mb)
        -1.0d0 * t1_A(ia, ba) * t2_A(ma, jb, na, mb) &
        +1.0d0 * t1_A(ja, ba) * t2_A(ma, ib, na, mb)
        if (ib /= nb .and. bb /= mb) then
        M2_A(ia,ja,na,ba) = M2_A(ia,ja,na,ba) &
        +1.0d0 * t1_A(ib, bb) * t2_A(ma, jb, na, mb)
        endif
        if (jb /= nb .and. bb /= mb) then
        M2_A(ia,ja,na,ba) = M2_A(ia,ja,na,ba) &
        -1.0d0 * t1_A(jb, bb) * t2_A(ma, ib, na, mb)
        endif
      enddo
    enddo
  enddo

  !! Deltas:((na, ba))
      do aa = i_aa, f_aa
!        if (aa == na) cycle
        ab = aa + cc_nVa
    do ja = i_ja, f_ja
!      if (ja == ma) cycle
      jb = ja + cc_nOa
  do ia = i_ia, f_ia
!    if (ia == ma) cycle
    ib = ia + cc_nOa
        M2_A(ia,ja,aa,na) = M2_A(ia,ja,aa,na) &
        !-1.0d0 * t1_A(ia, aa) * t2_A(ma, jb, na, mb) &
        !+1.0d0 * t1_A(ja, aa) * t2_A(ma, ib, na, mb) &
        !+1.0d0 * t1_A(ib, ab) * t2_A(ma, jb, na, mb) &
        !-1.0d0 * t1_A(jb, ab) * t2_A(ma, ib, na, mb)
        +1.0d0 * t1_A(ia, aa) * t2_A(ma, jb, na, mb) &
        -1.0d0 * t1_A(ja, aa) * t2_A(ma, ib, na, mb) 
        if (ib /= nb .and. ab /= mb) then
        M2_A(ia,ja,aa,na) = M2_A(ia,ja,aa,na) &
        -1.0d0 * t1_A(ib, ab) * t2_A(ma, jb, na, mb) 
        endif
        if (jb /= nb .and. ab /= mb ) then
        M2_A(ia,ja,aa,na) = M2_A(ia,ja,aa,na) &
        +1.0d0 * t1_A(jb, ab) * t2_A(ma, ib, na, mb)
        endif
      enddo
    enddo
  enddo

  ! ### Spin case: i_a, j_b, a_a, b_b ###
  !! Deltas:((na, aa))
      do bb = i_bb, f_bb
!        if (bb == mb) cycle
        ba = bb - cc_nVa
    do jb = i_jb, f_jb
!      if (jb == nb) cycle
      ja = jb - cc_nOa
  do ia = i_ia, f_ia
!    if (ia == ma) cycle
    ib = ia + cc_nOa
        M2_A(ia,jb,na,bb) = M2_A(ia,jb,na,bb) &
        !-1.0d0 * t1_A(jb, bb) * t2_A(ma, ib, na, mb) &
        !+1.0d0 * t1_A(ja, ba) * t2_A(ma, ib, na, mb)
        +1.0d0 * t1_A(jb, bb) * t2_A(ma, ib, na, mb) 
        if (ja /= ma .and. ba /= na) then
        M2_A(ia,jb,na,bb) = M2_A(ia,jb,na,bb) &
        -1.0d0 * t1_A(ja, ba) * t2_A(ma, ib, na, mb)
        endif
      enddo
    enddo
  enddo

  !! Deltas:((mb, bb))
      do aa = i_aa, f_aa
!        if (aa == na) cycle
        ab = aa + cc_nVa
    do jb = i_jb, f_jb
!      if (jb == nb) cycle
      ja = jb - cc_nOa
  do ia = i_ia, f_ia
!    if (ia == ma) cycle
    ib = ia + cc_nOa
        M2_A(ia,jb,aa,mb) = M2_A(ia,jb,aa,mb) &
        !-1.0d0 * t1_A(ia, aa) * t2_A(ja, nb, na, mb) &
        !+1.0d0 * t1_A(ib, ab) * t2_A(ja, nb, na, mb)
        +1.0d0 * t1_A(ia, aa) * t2_A(ja, nb, na, mb) 
        if (ib /= nb .and. ab /= mb) then
        M2_A(ia,jb,aa,mb) = M2_A(ia,jb,aa,mb) &
        -1.0d0 * t1_A(ib, ab) * t2_A(ja, nb, na, mb)
        endif
      enddo
    enddo
  enddo

  ! ### Spin case: i_a, j_b, a_b, b_a ###
  !! Deltas:((mb, ab))
      do ba = i_ba, f_ba
!        if (ba == na) cycle
        bb = ba + cc_nVa
    do jb = i_jb, f_jb
!      if (jb == nb) cycle
      ja = jb - cc_nOa
  do ia = i_ia, f_ia
!    if (ia == ma) cycle
    ib = ia + cc_nOa
        M2_A(ia,jb,mb,ba) = M2_A(ia,jb,mb,ba) &
        !+1.0d0 * t1_A(ia, ba) * t2_A(ja, nb, na, mb) &
        !-1.0d0 * t1_A(ib, bb) * t2_A(ja, nb, na, mb)
        -1.0d0 * t1_A(ia, ba) * t2_A(ja, nb, na, mb) 
        if (ib /= nb .and. bb /= mb) then
        M2_A(ia,jb,mb,ba) = M2_A(ia,jb,mb,ba) &
        +1.0d0 * t1_A(ib, bb) * t2_A(ja, nb, na, mb)
        endif
      enddo
    enddo
  enddo

  !! Deltas:((na, ba))
      do ab = i_ab, f_ab
!        if (ab == mb) cycle
        aa = ab - cc_nVa
    do jb = i_jb, f_jb
!      if (jb == nb) cycle
      ja = jb - cc_nOa
  do ia = i_ia, f_ia
!    if (ia == ma) cycle
    ib = ia + cc_nOa
        M2_A(ia,jb,ab,na) = M2_A(ia,jb,ab,na) &
        !+1.0d0 * t1_A(jb, ab) * t2_A(ma, ib, na, mb) &
        !-1.0d0 * t1_A(ja, aa) * t2_A(ma, ib, na, mb)
        -1.0d0 * t1_A(jb, ab) * t2_A(ma, ib, na, mb)
        if (ja /= ma .and. aa /= na) then
        M2_A(ia,jb,ab,na) = M2_A(ia,jb,ab,na) &
        +1.0d0 * t1_A(ja, aa) * t2_A(ma, ib, na, mb)
        endif
      enddo
    enddo
  enddo

  ! ### Spin case: i_b, j_a, a_a, b_b ###
  !! Deltas:((na, aa))
      do bb = i_bb, f_bb
!        if (bb == mb) cycle
        ba = bb - cc_nVa
    do ja = i_ja, f_ja
!      if (ja == ma) cycle
      jb = ja + cc_nOa
  do ib = i_ib, f_ib
!    if (ib == nb) cycle
    ia = ib - cc_nOa
        M2_A(ib,ja,na,bb) = M2_A(ib,ja,na,bb) &
        !+1.0d0 * t1_A(ib, bb) * t2_A(ma, jb, na, mb) &
        !-1.0d0 * t1_A(ia, ba) * t2_A(ma, jb, na, mb)
        -1.0d0 * t1_A(ib, bb) * t2_A(ma, jb, na, mb)
        if (ia /= ma .and. ba /= na) then
        M2_A(ib,ja,na,bb) = M2_A(ib,ja,na,bb) &
        +1.0d0 * t1_A(ia, ba) * t2_A(ma, jb, na, mb)
        endif
      enddo
    enddo
  enddo

  !! Deltas:((mb, bb))
      do aa = i_aa, f_aa
!        if (aa == na) cycle
        ab = aa + cc_nVa
    do ja = i_ja, f_ja
!      if (ja == ma) cycle
      jb = ja + cc_nOa
  do ib = i_ib, f_ib
!    if (ib == nb) cycle
    ia = ib - cc_nOa
        M2_A(ib,ja,aa,mb) = M2_A(ib,ja,aa,mb) &
        !+1.0d0 * t1_A(ja, aa) * t2_A(ia, nb, na, mb) &
        !-1.0d0 * t1_A(jb, ab) * t2_A(ia, nb, na, mb)
        -1.0d0 * t1_A(ja, aa) * t2_A(ia, nb, na, mb)
        if (jb /= nb .and. ab /= mb) then
        M2_A(ib,ja,aa,mb) = M2_A(ib,ja,aa,mb) &
        +1.0d0 * t1_A(jb, ab) * t2_A(ia, nb, na, mb)
        endif
      enddo
    enddo
  enddo

  ! ### Spin case: i_b, j_a, a_b, b_a ###
  !! Deltas:((mb, ab))
      do ba = i_ba, f_ba
!        if (ba == na) cycle
        bb = ba + cc_nVa
    do ja = i_ja, f_ja
!      if (ja == ma) cycle
      jb = ja + cc_nOa
  do ib = i_ib, f_ib
!    if (ib == nb) cycle
    ia = ib - cc_nOa
        M2_A(ib,ja,mb,ba) = M2_A(ib,ja,mb,ba) &
        !-1.0d0 * t1_A(ja, ba) * t2_A(ia, nb, na, mb) &
        !+1.0d0 * t1_A(jb, bb) * t2_A(ia, nb, na, mb)
        +1.0d0 * t1_A(ja, ba) * t2_A(ia, nb, na, mb)
        if (jb /= nb .and. bb /= mb) then
        M2_A(ib,ja,mb,ba) = M2_A(ib,ja,mb,ba) &
        -1.0d0 * t1_A(jb, bb) * t2_A(ia, nb, na, mb)
        endif
      enddo
    enddo
  enddo

  !! Deltas:((na, ba))
      do ab = i_ab, f_ab
!        if (ab == mb) cycle
        aa = ab - cc_nVa
    do ja = i_ja, f_ja
!      if (ja == ma) cycle
      jb = ja + cc_nOa
  do ib = i_ib, f_ib
!    if (ib == nb) cycle
    ia = ib - cc_nOa
        M2_A(ib,ja,ab,na) = M2_A(ib,ja,ab,na) &
        !-1.0d0 * t1_A(ib, ab) * t2_A(ma, jb, na, mb) &
        !+1.0d0 * t1_A(ia, aa) * t2_A(ma, jb, na, mb)
        +1.0d0 * t1_A(ib, ab) * t2_A(ma, jb, na, mb)
        if (ia /= ma .and. aa /= na) then
        M2_A(ib,ja,ab,na) = M2_A(ib,ja,ab,na) &
        -1.0d0 * t1_A(ia, aa) * t2_A(ma, jb, na, mb)
        endif
      enddo
    enddo
  enddo

  ! ### Spin case: i_b, j_b, a_b, b_b ###
  !! Deltas:((mb, ab))
      do bb = i_bb, f_bb
!        if (bb == mb) cycle 
        ba = bb - cc_nVa
    do jb = i_jb, f_jb
!      if (jb == nb) cycle 
      ja = jb - cc_nOa
  do ib = i_ib, f_ib
!    if (ib == nb) cycle 
    ia = ib - cc_nOa
        M2_A(ib,jb,mb,bb) = M2_A(ib,jb,mb,bb) & 
        !+1.0d0 * t1_A(ib, bb) * t2_A(ja, nb, na, mb) & 
        !-1.0d0 * t1_A(jb, bb) * t2_A(ia, nb, na, mb) & 
        !-1.0d0 * t1_A(ia, ba) * t2_A(ja, nb, na, mb) & 
        !+1.0d0 * t1_A(ja, ba) * t2_A(ia, nb, na, mb)
        -1.0d0 * t1_A(ib, bb) * t2_A(ja, nb, na, mb) & 
        +1.0d0 * t1_A(jb, bb) * t2_A(ia, nb, na, mb) 
        if (ia /= ma .and. ba /= na) then
        M2_A(ib,jb,mb,bb) = M2_A(ib,jb,mb,bb) & 
        +1.0d0 * t1_A(ia, ba) * t2_A(ja, nb, na, mb)  
        endif
        if (ja /= ma .and. ba /= na) then
        M2_A(ib,jb,mb,bb) = M2_A(ib,jb,mb,bb) & 
        -1.0d0 * t1_A(ja, ba) * t2_A(ia, nb, na, mb)
        endif
      enddo
    enddo
  enddo

  !! Deltas:((mb, bb))
      do ab = i_ab, f_ab
!        if (ab == mb) cycle 
        aa = ab - cc_nVa
    do jb = i_jb, f_jb
!      if (jb == nb) cycle 
      ja = jb - cc_nOa
  do ib = i_ib, f_ib
!    if (ib == nb) cycle 
    ia = ib - cc_nOa
        M2_A(ib,jb,ab,mb) = M2_A(ib,jb,ab,mb) & 
        !-1.0d0 * t1_A(ib, ab) * t2_A(ja, nb, na, mb) & 
        !+1.0d0 * t1_A(jb, ab) * t2_A(ia, nb, na, mb) & 
        !+1.0d0 * t1_A(ia, aa) * t2_A(ja, nb, na, mb) & 
        !-1.0d0 * t1_A(ja, aa) * t2_A(ia, nb, na, mb)
        +1.0d0 * t1_A(ib, ab) * t2_A(ja, nb, na, mb) & 
        -1.0d0 * t1_A(jb, ab) * t2_A(ia, nb, na, mb) 
        if (ia /= ma .and. aa /= na) then
        M2_A(ib,jb,ab,mb) = M2_A(ib,jb,ab,mb) & 
        -1.0d0 * t1_A(ia, aa) * t2_A(ja, nb, na, mb)  
        endif
        if (ja /= ma .and. aa /= na) then
        M2_A(ib,jb,ab,mb) = M2_A(ib,jb,ab,mb) & 
        +1.0d0 * t1_A(ja, aa) * t2_A(ia, nb, na, mb)
        endif
      enddo
    enddo
  enddo

end

