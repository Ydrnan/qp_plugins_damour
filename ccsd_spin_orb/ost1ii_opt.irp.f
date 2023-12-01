subroutine ost1ii_opt(nO,nV,det,t1_A,t2_A,M1_A,M2_A)

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

    do aa = i_aa, f_aa
      if (aa == na) cycle
      ab = aa + cc_nVa
  do ia = i_ia, f_ia
    if (ia == ma) cycle
    ib = ia + cc_nOa
      M1_A(ia,aa) = M1_A(ia,aa) &
      -1.0d0 * t1_A(nb, ab) * t2_A(ma, ib, na, mb) &
      -1.0d0 * t1_A(ib, mb) * t2_A(ma, nb, na, ab)
    enddo
  enddo

  ! ### Spin case: i_b, a_b ###

    do ab = i_ab, f_ab
      if (ab == mb) cycle
      aa = ab - cc_nVa
  do ib = i_ib, f_ib
    if (ib == nb) cycle
    ia = ib - cc_nOa
      M1_A(ib,ab) = M1_A(ib,ab) &
      -1.0d0 * t1_A(ma, aa) * t2_A(ia, nb, na, mb) &
      -1.0d0 * t1_A(ia, na) * t2_A(ma, nb, aa, mb)
    enddo
  enddo

!    ! ### Spin case: i_a, j_a, a_a, b_a ###
!
!        do ba = i_ba, f_ba
!!          if (ba == na) cycle
!          bb = ba + cc_nVa
!      do aa = i_aa, f_aa
!!        if (aa == na) cycle
!        ab = aa + cc_nVa
!    do ja = i_ja, f_ja
!!      if (ja == ma) cycle
!      jb = ja + cc_nOa
!  do ia = i_ia, f_ia
!!    if (ia == ma) cycle
!    ib = ia + cc_nOa
!          M2_A(ia,ja,aa,ba) = M2_A(ia,ja,aa,ba) &
!          !-1.0d0 * t1_A(ia, aa) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) &
!          !-1.0d0 * t1_A(ia, aa) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) &
!          !+1.0d0 * t1_A(ja, aa) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) &
!          !+1.0d0 * t1_A(ja, aa) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) &
!          !+1.0d0 * t1_A(ia, ba) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) &
!          !+1.0d0 * t1_A(ia, ba) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) &
!          !-1.0d0 * t1_A(ja, ba) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) &
!          !-1.0d0 * t1_A(ja, ba) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) &
!          !+1.0d0 * t1_A(ib, ab) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) &
!          !+1.0d0 * t1_A(ib, ab) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) &
!          !-1.0d0 * t1_A(jb, ab) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) &
!          !-1.0d0 * t1_A(jb, ab) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) &
!          !-1.0d0 * t1_A(ib, bb) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) &
!          !-1.0d0 * t1_A(ib, bb) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) &
!          !+1.0d0 * t1_A(jb, bb) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) &
!          !+1.0d0 * t1_A(jb, bb) * t1_A(ib, mb) * t2_A(ma, nb, na, ab)
!        
!          +1.0d0 * t1_A(ia, aa) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) &
!          +1.0d0 * t1_A(ia, aa) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) &
!          -1.0d0 * t1_A(ja, aa) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) &
!          -1.0d0 * t1_A(ja, aa) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) &
!          -1.0d0 * t1_A(ia, ba) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) &
!          -1.0d0 * t1_A(ia, ba) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) &
!          +1.0d0 * t1_A(ja, ba) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) &
!          +1.0d0 * t1_A(ja, ba) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) 
!          if (ib /= nb .and. ab /= mb) then
!          ! .and. jb /= nb .and. bb /= mb) then
!          M2_A(ia,ja,aa,ba) = M2_A(ia,ja,aa,ba) &
!          -1.0d0 * t1_A(ib, ab) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) &
!          -1.0d0 * t1_A(ib, ab) * t1_A(jb, mb) * t2_A(ma, nb, na, bb)
!          endif
!          if (jb /= nb .and. ab /= mb) then
!          M2_A(ia,ja,aa,ba) = M2_A(ia,ja,aa,ba) &
!          +1.0d0 * t1_A(jb, ab) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) &
!          +1.0d0 * t1_A(jb, ab) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) 
!          endif
!          if (ib /= nb .and. bb /= mb) then
!          M2_A(ia,ja,aa,ba) = M2_A(ia,ja,aa,ba) &
!          +1.0d0 * t1_A(ib, bb) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) &
!          +1.0d0 * t1_A(ib, bb) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) 
!          endif
!          if (jb /= nb .and. bb /= mb) then
!          M2_A(ia,ja,aa,ba) = M2_A(ia,ja,aa,ba) &
!          -1.0d0 * t1_A(jb, bb) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) &
!          -1.0d0 * t1_A(jb, bb) * t1_A(ib, mb) * t2_A(ma, nb, na, ab)
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!
!!  do ba = i_ba, f_ba
!!    bb = ba + cc_nVa
!!    do aa = i_aa, f_aa
!!      ab = aa + cc_nVa
!!      do ja = i_ja, f_ja
!!        jb = ja + cc_nOa
!!        do ia = i_ia, f_ia
!!          ib = ia + cc_nOa
!!          M2_A(ia,ja,aa,ba) = M2_A(ia,ja,aa,ba) &
!!            + (t2_A(ib,ma,mb,na) * t1_A(nb,ab) + t2_A(nb,ma,ab,na) * t1_A(ib,mb)) * (t1_A(ja,ba)) &
!!            - (t2_A(jb,ma,mb,na) * t1_A(nb,ab) + t2_A(nb,ma,ab,na) * t1_A(jb,mb)) * (t1_A(ia,ba)) &
!!            - (t2_A(ib,ma,mb,na) * t1_A(nb,bb) + t2_A(nb,ma,bb,na) * t1_A(ib,mb)) * (t1_A(ja,aa)) &
!!            + (t2_A(jb,ma,mb,na) * t1_A(nb,bb) + t2_A(nb,ma,bb,na) * t1_A(jb,mb)) * (t1_A(ia,aa))
!!          if (ib /= nb .and. jb /= nb .and. ab /= mb .and. bb /= mb) then
!!          M2_A(ia,ja,aa,ba) = M2_A(ia,ja,aa,ba) &
!!            + (t2_A(ib,ma,mb,na) * t1_A(nb,ab) + t2_A(nb,ma,ab,na) * t1_A(ib,mb)) * (- t1_A(jb,bb)) &
!!            - (t2_A(jb,ma,mb,na) * t1_A(nb,ab) + t2_A(nb,ma,ab,na) * t1_A(jb,mb)) * (- t1_A(ib,bb)) &
!!            - (t2_A(ib,ma,mb,na) * t1_A(nb,bb) + t2_A(nb,ma,bb,na) * t1_A(ib,mb)) * (- t1_A(jb,ab)) &
!!            + (t2_A(jb,ma,mb,na) * t1_A(nb,bb) + t2_A(nb,ma,bb,na) * t1_A(jb,mb)) * (- t1_A(ib,ab))
!!          endif
!!        enddo
!!      enddo
!!    enddo
!!  enddo
!
!
!  ! ### Spin case: i_a, j_b, a_a, b_b ###
!
!        do bb = i_bb, f_bb
!!          if (bb == mb) cycle
!          ba = bb - cc_nVa
!      do aa = i_aa, f_aa
!!        if (aa == na) cycle
!        ab = aa + cc_nVa
!    do jb = i_jb, f_jb
!!      if (jb == nb) cycle
!      ja = jb - cc_nOa
!  do ia = i_ia, f_ia
!!    if (ia == ma) cycle
!    ib = ia + cc_nOa
!          M2_A(ia,jb,aa,bb) = M2_A(ia,jb,aa,bb) &
!          !-1.0d0 * t1_A(ia, aa) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) &
!          !-1.0d0 * t1_A(ia, aa) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) &
!          !-1.0d0 * t1_A(jb, bb) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) &
!          !-1.0d0 * t1_A(jb, bb) * t1_A(ib, mb) * t2_A(ma, nb, na, ab) &
!          !+1.0d0 * t1_A(ib, ab) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) &
!          !+1.0d0 * t1_A(ib, ab) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) &
!          !+1.0d0 * t1_A(ja, ba) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) &
!          !+1.0d0 * t1_A(ja, ba) * t1_A(ib, mb) * t2_A(ma, nb, na, ab)
!
!          +1.0d0 * t1_A(ia, aa) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) &
!          +1.0d0 * t1_A(ia, aa) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) &
!          +1.0d0 * t1_A(jb, bb) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) &
!          +1.0d0 * t1_A(jb, bb) * t1_A(ib, mb) * t2_A(ma, nb, na, ab)
!          if (ib /= nb .and. ab /= mb .and. ja /= ma .and. ba /= na) then
!          M2_A(ia,jb,aa,bb) = M2_A(ia,jb,aa,bb) &
!          -1.0d0 * t1_A(ib, ab) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) &
!          -1.0d0 * t1_A(ib, ab) * t1_A(ja, na) * t2_A(ma, nb, ba, mb)
!          !endif
!          !if (ja /= ma .and. ba /= na) then
!          M2_A(ia,jb,aa,bb) = M2_A(ia,jb,aa,bb) &
!          -1.0d0 * t1_A(ja, ba) * t1_A(nb, ab) * t2_A(ma, ib, na, mb) &
!          -1.0d0 * t1_A(ja, ba) * t1_A(ib, mb) * t2_A(ma, nb, na, ab)
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!
!    ! ### Spin case: i_a, j_b, a_b, b_a ###
!
!        do ba = i_ba, f_ba
!!          if (ba == na) cycle
!          bb = ba + cc_nVa
!      do ab = i_ab, f_ab
!!        if (ab == mb) cycle
!        aa = ab - cc_nVa
!    do jb = i_jb, f_jb
!!      if (jb == nb) cycle
!      ja = jb - cc_nOa
!  do ia = i_ia, f_ia
!!    if (ia == ma) cycle
!    ib = ia + cc_nOa
!          M2_A(ia,jb,ab,ba) = M2_A(ia,jb,ab,ba) &
!          !+1.0d0 * t1_A(jb, ab) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) &
!          !+1.0d0 * t1_A(jb, ab) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) &
!          !+1.0d0 * t1_A(ia, ba) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) &
!          !+1.0d0 * t1_A(ia, ba) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) &
!          !-1.0d0 * t1_A(ja, aa) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) &
!          !-1.0d0 * t1_A(ja, aa) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) &
!          !-1.0d0 * t1_A(ib, bb) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) &
!          !-1.0d0 * t1_A(ib, bb) * t1_A(ja, na) * t2_A(ma, nb, aa, mb)
!
!          -1.0d0 * t1_A(jb, ab) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) &
!          -1.0d0 * t1_A(jb, ab) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) &
!          -1.0d0 * t1_A(ia, ba) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) &
!          -1.0d0 * t1_A(ia, ba) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) 
!          if (ja /= ma .and. aa /= na .and. ib /= nb .and. bb /= mb) then
!          M2_A(ia,jb,ab,ba) = M2_A(ia,jb,ab,ba) &
!          +1.0d0 * t1_A(ja, aa) * t1_A(nb, bb) * t2_A(ma, ib, na, mb) &
!          +1.0d0 * t1_A(ja, aa) * t1_A(ib, mb) * t2_A(ma, nb, na, bb) 
!          !endif
!          !if (ib /= nb .and. bb /= mb) then
!          M2_A(ia,jb,ab,ba) = M2_A(ia,jb,ab,ba) &
!          +1.0d0 * t1_A(ib, bb) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) &
!          +1.0d0 * t1_A(ib, bb) * t1_A(ja, na) * t2_A(ma, nb, aa, mb)
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!
! ! ### Spin case: i_b, j_a, a_a, b_b ###
!
!        do bb = i_bb, f_bb
!!          if (bb == mb) cycle
!          ba = bb - cc_nVa
!      do aa = i_aa, f_aa
!!        if (aa == na) cycle
!        ab = aa + cc_nVa
!    do ja = i_ja, f_ja
!!      if (ja == ma) cycle
!      jb = ja + cc_nOa
!  do ib = i_ib, f_ib
!!    if (ib == nb) cycle
!    ia = ib - cc_nOa
!          M2_A(ib,ja,aa,bb) = M2_A(ib,ja,aa,bb) &
!          !+1.0d0 * t1_A(ja, aa) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) &
!          !+1.0d0 * t1_A(ja, aa) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) &
!          !+1.0d0 * t1_A(ib, bb) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) &
!          !+1.0d0 * t1_A(ib, bb) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) &
!          !-1.0d0 * t1_A(jb, ab) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) &
!          !-1.0d0 * t1_A(jb, ab) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) &
!          !-1.0d0 * t1_A(ia, ba) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) &
!          !-1.0d0 * t1_A(ia, ba) * t1_A(jb, mb) * t2_A(ma, nb, na, ab)
!
!          -1.0d0 * t1_A(ja, aa) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) &
!          -1.0d0 * t1_A(ja, aa) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) &
!          -1.0d0 * t1_A(ib, bb) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) &
!          -1.0d0 * t1_A(ib, bb) * t1_A(jb, mb) * t2_A(ma, nb, na, ab) 
!          if (jb /= nb .and. ab /= mb .and. ia /= ma .and. ba /= na) then
!          M2_A(ib,ja,aa,bb) = M2_A(ib,ja,aa,bb) &
!          +1.0d0 * t1_A(jb, ab) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) &
!          +1.0d0 * t1_A(jb, ab) * t1_A(ia, na) * t2_A(ma, nb, ba, mb)
!          !endif
!          !if (ia /= ma .and. ba /= na) then
!          M2_A(ib,ja,aa,bb) = M2_A(ib,ja,aa,bb) &
!          +1.0d0 * t1_A(ia, ba) * t1_A(nb, ab) * t2_A(ma, jb, na, mb) &
!          +1.0d0 * t1_A(ia, ba) * t1_A(jb, mb) * t2_A(ma, nb, na, ab)
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!
!  ! ### Spin case: i_b, j_a, a_b, b_a ###
!
!        do ba = i_ba, f_ba
!!          if (ba == na) cycle
!          bb = ba + cc_nVa
!      do ab = i_ab, f_ab
!!        if (ab == mb) cycle
!        aa = ab - cc_nVa
!    do ja = i_ja, f_ja
!!      if (ja == ma) cycle
!      jb = ja + cc_nOa
!  do ib = i_ib, f_ib
!!    if (ib == nb) cycle
!    ia = ib - cc_nOa
!          M2_A(ib,ja,ab,ba) = M2_A(ib,ja,ab,ba) &
!          !-1.0d0 * t1_A(ib, ab) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) &
!          !-1.0d0 * t1_A(ib, ab) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) &
!          !-1.0d0 * t1_A(ja, ba) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) &
!          !-1.0d0 * t1_A(ja, ba) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) &
!          !+1.0d0 * t1_A(ia, aa) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) &
!          !+1.0d0 * t1_A(ia, aa) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) &
!          !+1.0d0 * t1_A(jb, bb) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) &
!          !+1.0d0 * t1_A(jb, bb) * t1_A(ia, na) * t2_A(ma, nb, aa, mb)
!
!          +1.0d0 * t1_A(ib, ab) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) &
!          +1.0d0 * t1_A(ib, ab) * t1_A(jb, mb) * t2_A(ma, nb, na, bb) &
!          +1.0d0 * t1_A(ja, ba) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) &
!          +1.0d0 * t1_A(ja, ba) * t1_A(ia, na) * t2_A(ma, nb, aa, mb)
!          if (ia /= ma .and. aa /= na .and. jb /= nb .and. bb /= mb) then
!          M2_A(ib,ja,ab,ba) = M2_A(ib,ja,ab,ba) &
!          -1.0d0 * t1_A(ia, aa) * t1_A(nb, bb) * t2_A(ma, jb, na, mb) &
!          -1.0d0 * t1_A(ia, aa) * t1_A(jb, mb) * t2_A(ma, nb, na, bb)
!          !endif
!          !if (jb /= nb .and. bb /= mb) then
!          M2_A(ib,ja,ab,ba) = M2_A(ib,ja,ab,ba) &
!          -1.0d0 * t1_A(jb, bb) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) &
!          -1.0d0 * t1_A(jb, bb) * t1_A(ia, na) * t2_A(ma, nb, aa, mb)
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!
!  ! ### Spin case: i_b, j_b, a_b, b_b ###
!
!        do bb = i_bb, f_bb
!!          if (bb == mb) cycle
!          ba = bb - cc_nVa
!      do ab = i_ab, f_ab
!!        if (ab == mb) cycle
!        aa = ab - cc_nVa
!    do jb = i_jb, f_jb
!!      if (jb == nb) cycle
!      ja = jb - cc_nOa
!  do ib = i_ib, f_ib
!!    if (ib == nb) cycle
!    ia = ib - cc_nOa
!          M2_A(ib,jb,ab,bb) = M2_A(ib,jb,ab,bb) &
!          !-1.0d0 * t1_A(ib, ab) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) &
!          !-1.0d0 * t1_A(ib, ab) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) &
!          !+1.0d0 * t1_A(jb, ab) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) &
!          !+1.0d0 * t1_A(jb, ab) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) &
!          !+1.0d0 * t1_A(ib, bb) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) &
!          !+1.0d0 * t1_A(ib, bb) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) &
!          !-1.0d0 * t1_A(jb, bb) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) &
!          !-1.0d0 * t1_A(jb, bb) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) &
!          !+1.0d0 * t1_A(ia, aa) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) &
!          !+1.0d0 * t1_A(ia, aa) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) &
!          !-1.0d0 * t1_A(ja, aa) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) &
!          !-1.0d0 * t1_A(ja, aa) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) &
!          !-1.0d0 * t1_A(ia, ba) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) &
!          !-1.0d0 * t1_A(ia, ba) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) &
!          !+1.0d0 * t1_A(ja, ba) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) &
!          !+1.0d0 * t1_A(ja, ba) * t1_A(ia, na) * t2_A(ma, nb, aa, mb)
!
!          +1.0d0 * t1_A(ib, ab) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) &
!          +1.0d0 * t1_A(ib, ab) * t1_A(ja, na) * t2_A(ma, nb, ba, mb) &
!          -1.0d0 * t1_A(jb, ab) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) &
!          -1.0d0 * t1_A(jb, ab) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) &
!          -1.0d0 * t1_A(ib, bb) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) &
!          -1.0d0 * t1_A(ib, bb) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) &
!          +1.0d0 * t1_A(jb, bb) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) &
!          +1.0d0 * t1_A(jb, bb) * t1_A(ia, na) * t2_A(ma, nb, aa, mb) 
!          if (ia /= ma .and. aa /= na) then 
!          !.and. ja /= ma .and. ba /= na) then
!          M2_A(ib,jb,ab,bb) = M2_A(ib,jb,ab,bb) &
!          -1.0d0 * t1_A(ia, aa) * t1_A(ma, ba) * t2_A(ja, nb, na, mb) &
!          -1.0d0 * t1_A(ia, aa) * t1_A(ja, na) * t2_A(ma, nb, ba, mb)
!          endif 
!          if (ja /= ma .and. aa /= na) then
!          M2_A(ib,jb,ab,bb) = M2_A(ib,jb,ab,bb) &
!          +1.0d0 * t1_A(ja, aa) * t1_A(ma, ba) * t2_A(ia, nb, na, mb) &
!          +1.0d0 * t1_A(ja, aa) * t1_A(ia, na) * t2_A(ma, nb, ba, mb) 
!          endif 
!          if (ia /= ma .and. ba /= na) then
!          M2_A(ib,jb,ab,bb) = M2_A(ib,jb,ab,bb) &
!          +1.0d0 * t1_A(ia, ba) * t1_A(ma, aa) * t2_A(ja, nb, na, mb) &
!          +1.0d0 * t1_A(ia, ba) * t1_A(ja, na) * t2_A(ma, nb, aa, mb) 
!          endif 
!          if (ja /= ma .and. ba /= na) then
!          M2_A(ib,jb,ab,bb) = M2_A(ib,jb,ab,bb) &
!          -1.0d0 * t1_A(ja, ba) * t1_A(ma, aa) * t2_A(ia, nb, na, mb) &
!          -1.0d0 * t1_A(ja, ba) * t1_A(ia, na) * t2_A(ma, nb, aa, mb)
!          endif 
!        enddo
!      enddo
!    enddo
!  enddo
!
!!  do ba = i_ba, f_ba
!!    bb = ba + cc_nVa
!!    do aa = i_aa, f_aa
!!      ab = aa + cc_nVa
!!      do ja = i_ja, f_ja
!!        jb = ja + cc_nOa
!!        do ia = i_ia, f_ia
!!          ib = ia + cc_nOa
!!          M2_A(ib,jb,ab,bb) = M2_A(ib,jb,ab,bb) &
!!            - (t2_A(ia,nb,na,mb) * t1_A(ma,aa) + t2_A(nb,ma,mb,aa) * t1_A(ia,na)) * (t1_A(jb,bb)) &
!!            + (t2_A(ja,nb,na,mb) * t1_A(ma,aa) + t2_A(nb,ma,mb,aa) * t1_A(ja,na)) * (t1_A(ib,bb)) &
!!            + (t2_A(ia,nb,na,mb) * t1_A(ma,ba) + t2_A(nb,ma,mb,ba) * t1_A(ia,na)) * (t1_A(jb,ab)) &
!!            - (t2_A(ja,nb,na,mb) * t1_A(ma,ba) + t2_A(nb,ma,mb,ba) * t1_A(ja,na)) * (t1_A(ib,ab)) 
!!            if (ja /= ma .and. ba /= na) then
!!          M2_A(ib,jb,ab,bb) = M2_A(ib,jb,ab,bb) &
!!            - (t2_A(ia,nb,na,mb) * t1_A(ma,aa) + t2_A(nb,ma,mb,aa) * t1_A(ia,na)) * (- t1_A(ja,ba)) 
!!            endif
!!            if (ia /= ma .and. ba /= na) then
!!          M2_A(ib,jb,ab,bb) = M2_A(ib,jb,ab,bb) &
!!            + (t2_A(ja,nb,na,mb) * t1_A(ma,aa) + t2_A(nb,ma,mb,aa) * t1_A(ja,na)) * (- t1_A(ia,ba))
!!            endif
!!            if (ja /= ma .and. aa /= na) then
!!          M2_A(ib,jb,ab,bb) = M2_A(ib,jb,ab,bb) &
!!            + (t2_A(ia,nb,na,mb) * t1_A(ma,ba) + t2_A(nb,ma,mb,ba) * t1_A(ia,na)) * (- t1_A(ja,aa))
!!            endif
!!            if (ia /= ma .and. aa /= na) then
!!          M2_A(ib,jb,ab,bb) = M2_A(ib,jb,ab,bb) &
!!            - (t2_A(ja,nb,na,mb) * t1_A(ma,ba) + t2_A(nb,ma,mb,ba) * t1_A(ja,na)) * (- t1_A(ia,aa))
!!            endif
!!        enddo
!!      enddo
!!    enddo
!!  enddo

end

