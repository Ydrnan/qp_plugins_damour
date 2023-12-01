subroutine ost1ia_opt(nO,nV,det,t1_A,t2_A,M1_A,M2_A)

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

  !! Deltas:((ma, ia))
  do aa = i_aa, f_aa
    if (aa == na) cycle
    ab = aa + cc_nVa
    M1_A(ma,aa) = M1_A(ma,aa) &
    +1.0d0 * t2_A(ma, nb, na, ab)
  enddo

  ! ### Spin case: i_b, a_b ###

  !! Deltas:((nb, ib))
  do ab = i_ab, f_ab
    if (ab == mb) cycle
    aa = ab - cc_nVa
    M1_A(nb,ab) = M1_A(nb,ab) &
    +1.0d0 * t2_A(ma, nb, aa, mb)
  enddo

!  ! ### Spin case: i_a, j_a, a_a, b_a ###
!
!  !! Deltas:((ma, ja))
!  do ba = i_ba, f_ba
!    bb = ba + cc_nVa
!    do aa = i_aa, f_aa
!      ab = aa + cc_nVa
!      do ia = i_ia, f_ia
!        ib = ia + cc_nOa
!     !   M2_A(ia,ma,aa,ba) = M2_A(ia,ma,aa,ba) &
!     !   !+1.0d0 * t1_A(ia, aa) * t2_A(ma, nb, na, bb) &
!     !   !-1.0d0 * t1_A(ia, ba) * t2_A(ma, nb, na, ab) &
!     !   !-1.0d0 * t1_A(ib, ab) * t2_A(ma, nb, na, bb) &
!     !   !+1.0d0 * t1_A(ib, bb) * t2_A(ma, nb, na, ab)
!
!     !   -1.0d0 * t1_A(ia, aa) * t2_A(ma, nb, na, bb) &
!     !   +1.0d0 * t1_A(ia, ba) * t2_A(ma, nb, na, ab) 
!        
!        if (ib /= nb .and. ab /= mb) then
!        M2_A(ia,ma,aa,ba) = M2_A(ia,ma,aa,ba) &
!        +1.0d0 * t1_A(ib, ab) * t2_A(ma, nb, na, bb)
!        endif
!        if (ib /= nb .and. bb /= mb) then
!        M2_A(ia,ma,aa,ba) = M2_A(ia,ma,aa,ba) &
!        -1.0d0 * t1_A(ib, bb) * t2_A(ma, nb, na, ab)
!        endif
!
!      enddo
!    enddo
!  enddo
!
!  !! Deltas:((ma, ia))
!      do ba = i_ba, f_ba
!        bb = ba + cc_nVa
!    do aa = i_aa, f_aa
!      ab = aa + cc_nVa
!  do ja = i_ja, f_ja
!    jb = ja + cc_nOa
!     !   M2_A(ma,ja,aa,ba) = M2_A(ma,ja,aa,ba) &
!     !   !-1.0d0 * t1_A(ja, aa) * t2_A(ma, nb, na, bb) &
!     !   !+1.0d0 * t1_A(ja, ba) * t2_A(ma, nb, na, ab) &
!     !   !+1.0d0 * t1_A(jb, ab) * t2_A(ma, nb, na, bb) &
!     !   !-1.0d0 * t1_A(jb, bb) * t2_A(ma, nb, na, ab)
!     !   +1.0d0 * t1_A(ja, aa) * t2_A(ma, nb, na, bb) &
!     !   -1.0d0 * t1_A(ja, ba) * t2_A(ma, nb, na, ab)
!
!        if (jb /= nb .and. ab /= mb) then
!        M2_A(ma,ja,aa,ba) = M2_A(ma,ja,aa,ba) &
!        -1.0d0 * t1_A(jb, ab) * t2_A(ma, nb, na, bb) 
!        endif
!        if (jb /= nb .and. bb /= mb) then
!        M2_A(ma,ja,aa,ba) = M2_A(ma,ja,aa,ba) &
!        +1.0d0 * t1_A(jb, bb) * t2_A(ma, nb, na, ab)
!        endif
!      enddo
!    enddo
!  enddo
!
!  ! ### Spin case: i_a, j_b, a_a, b_b ###
!
!  !! Deltas:((nb, jb))
!      do bb = i_bb, f_bb
!        ba = bb - cc_nVa
!    do aa = i_aa, f_aa
!      ab = aa + cc_nVa
!  do ia = i_ia, f_ia
!    ib = ia + cc_nOa
!      !  M2_A(ia,nb,aa,bb) = M2_A(ia,nb,aa,bb) &
!      !  !+1.0d0 * t1_A(ia, aa) * t2_A(ma, nb, ba, mb) &
!      !  !-1.0d0 * t1_A(ib, ab) * t2_A(ma, nb, ba, mb)
!      !  -1.0d0 * t1_A(ia, aa) * t2_A(ma, nb, ba, mb) 
!        if (ib /= nb .and. ab /= mb) then
!        M2_A(ia,nb,aa,bb) = M2_A(ia,nb,aa,bb) &
!        +1.0d0 * t1_A(ib, ab) * t2_A(ma, nb, ba, mb)
!        endif
!      enddo
!    enddo
!  enddo
!
!  !! Deltas:((ma, ia))
!      do bb = i_bb, f_bb
!        ba = bb - cc_nVa
!    do aa = i_aa, f_aa
!      ab = aa + cc_nVa
!  do jb = i_jb, f_jb
!    ja = jb - cc_nOa
!      !  M2_A(ma,jb,aa,bb) = M2_A(ma,jb,aa,bb) &
!      !  !+1.0d0 * t1_A(jb, bb) * t2_A(ma, nb, na, ab) &
!      !  !-1.0d0 * t1_A(ja, ba) * t2_A(ma, nb, na, ab)
!      !  -1.0d0 * t1_A(jb, bb) * t2_A(ma, nb, na, ab) 
!        if (ja /= ma .and. ba /= na) then
!        M2_A(ma,jb,aa,bb) = M2_A(ma,jb,aa,bb) &
!        +1.0d0 * t1_A(ja, ba) * t2_A(ma, nb, na, ab)
!        endif
!      enddo
!    enddo
!  enddo
!
!  ! ### Spin case: i_a, j_b, a_b, b_a ###
!  !! Deltas:((nb, jb))
!      do ba = i_ba, f_ba
!        bb = ba + cc_nVa
!    do ab = i_ab, f_ab
!      aa = ab - cc_nVa
!  do ia = i_ia, f_ia
!    ib = ia + cc_nOa
!      !  M2_A(ia,nb,ab,ba) = M2_A(ia,nb,ab,ba) &
!      !  !-1.0d0 * t1_A(ia, ba) * t2_A(ma, nb, aa, mb) &
!      !  !+1.0d0 * t1_A(ib, bb) * t2_A(ma, nb, aa, mb)
!      !  +1.0d0 * t1_A(ia, ba) * t2_A(ma, nb, aa, mb) 
!        if (ib /= nb .and. bb /= mb) then
!        M2_A(ia,nb,ab,ba) = M2_A(ia,nb,ab,ba) &
!        -1.0d0 * t1_A(ib, bb) * t2_A(ma, nb, aa, mb)
!        endif
!      enddo
!    enddo
!  enddo
!
!  !! Deltas:((ma, ia))
!      do ba = i_ba, f_ba
!        bb = ba + cc_nVa
!    do ab = i_ab, f_ab
!      aa = ab - cc_nVa
!  do jb = i_jb, f_jb
!    ja = jb - cc_nOa
!      !  M2_A(ma,jb,ab,ba) = M2_A(ma,jb,ab,ba) &
!      !  !-1.0d0 * t1_A(jb, ab) * t2_A(ma, nb, na, bb) &
!      !  !+1.0d0 * t1_A(ja, aa) * t2_A(ma, nb, na, bb)
!      !  +1.0d0 * t1_A(jb, ab) * t2_A(ma, nb, na, bb) 
!        if (ja /= ma .and. aa /= na) then
!        M2_A(ma,jb,ab,ba) = M2_A(ma,jb,ab,ba) &
!        -1.0d0 * t1_A(ja, aa) * t2_A(ma, nb, na, bb)
!        endif
!      enddo
!    enddo
!  enddo
!
!  ! ### Spin case: i_b, j_a, a_a, b_b ###
!  !! Deltas:((ma, ja))
!      do bb = i_bb, f_bb
!        ba = bb - cc_nVa
!    do aa = i_aa, f_aa
!      ab = aa + cc_nVa
!  do ib = i_ib, f_ib
!    ia = ib - cc_nOa
!      !  M2_A(ib,ma,aa,bb) = M2_A(ib,ma,aa,bb) &
!      !  !-1.0d0 * t1_A(ib, bb) * t2_A(ma, nb, na, ab) &
!      !  !+1.0d0 * t1_A(ia, ba) * t2_A(ma, nb, na, ab)
!      !  +1.0d0 * t1_A(ib, bb) * t2_A(ma, nb, na, ab) 
!        if (ia /= ma .and. ba /= na) then
!        M2_A(ib,ma,aa,bb) = M2_A(ib,ma,aa,bb) &
!        -1.0d0 * t1_A(ia, ba) * t2_A(ma, nb, na, ab)
!        endif
!      enddo
!    enddo
!  enddo
!
!  !! Deltas:((nb, ib))
!      do bb = i_bb, f_bb
!        ba = bb - cc_nVa
!    do aa = i_aa, f_aa
!      ab = aa + cc_nVa
!  do ja = i_ja, f_ja
!    jb = ja + cc_nOa
!      !  M2_A(nb,ja,aa,bb) = M2_A(nb,ja,aa,bb) &
!      !  !-1.0d0 * t1_A(ja, aa) * t2_A(ma, nb, ba, mb) &
!      !  !+1.0d0 * t1_A(jb, ab) * t2_A(ma, nb, ba, mb)
!      !  +1.0d0 * t1_A(ja, aa) * t2_A(ma, nb, ba, mb)
!        if (jb /= nb .and. ab /= mb) then
!        M2_A(nb,ja,aa,bb) = M2_A(nb,ja,aa,bb) &
!        -1.0d0 * t1_A(jb, ab) * t2_A(ma, nb, ba, mb)
!        endif
!      enddo
!    enddo
!  enddo
!
!  ! ### Spin case: i_b, j_a, a_b, b_a ###
!  !! Deltas:((ma, ja))
!      do ba = i_ba, f_ba
!        bb = ba + cc_nVa
!    do ab = i_ab, f_ab
!      aa = ab - cc_nVa
!  do ib = i_ib, f_ib
!    ia = ib - cc_nOa
!      !  M2_A(ib,ma,ab,ba) = M2_A(ib,ma,ab,ba) &
!      !  !+1.0d0 * t1_A(ib, ab) * t2_A(ma, nb, na, bb) &
!      !  !-1.0d0 * t1_A(ia, aa) * t2_A(ma, nb, na, bb)
!      !  -1.0d0 * t1_A(ib, ab) * t2_A(ma, nb, na, bb)
!        if (ia /= ma .and. aa /= na) then
!        M2_A(ib,ma,ab,ba) = M2_A(ib,ma,ab,ba) &
!        +1.0d0 * t1_A(ia, aa) * t2_A(ma, nb, na, bb)
!        endif
!      enddo
!    enddo
!  enddo
!
!  !! Deltas:((nb, ib))
!      do ba = i_ba, f_ba
!        bb = ba + cc_nVa
!    do ab = i_ab, f_ab
!      aa = ab - cc_nVa
!  do ja = i_ja, f_ja
!    jb = ja + cc_nOa
!      !  M2_A(nb,ja,ab,ba) = M2_A(nb,ja,ab,ba) &
!      !  !+1.0d0 * t1_A(ja, ba) * t2_A(ma, nb, aa, mb) &
!      !  !-1.0d0 * t1_A(jb, bb) * t2_A(ma, nb, aa, mb)
!      !  -1.0d0 * t1_A(ja, ba) * t2_A(ma, nb, aa, mb) 
!        if (jb /= nb .and. bb /= mb) then
!        M2_A(nb,ja,ab,ba) = M2_A(nb,ja,ab,ba) &
!        +1.0d0 * t1_A(jb, bb) * t2_A(ma, nb, aa, mb)
!        endif
!      enddo
!    enddo
!  enddo
!
!  ! ### Spin case: i_b, j_b, a_b, b_b ###
!
!  !! Deltas:((nb, jb))
!      do bb = i_bb, f_bb
!        ba = bb - cc_nVa
!    do ab = i_ab, f_ab
!      aa = ab - cc_nVa
!  do ib = i_ib, f_ib
!    ia = ib - cc_nOa
!     !   M2_A(ib,nb,ab,bb) = M2_A(ib,nb,ab,bb) & 
!     !   !+1.0d0 * t1_A(ib, ab) * t2_A(ma, nb, ba, mb) & 
!     !   !-1.0d0 * t1_A(ib, bb) * t2_A(ma, nb, aa, mb) & 
!     !   !-1.0d0 * t1_A(ia, aa) * t2_A(ma, nb, ba, mb) & 
!     !   !+1.0d0 * t1_A(ia, ba) * t2_A(ma, nb, aa, mb)
!     !   -1.0d0 * t1_A(ib, ab) * t2_A(ma, nb, ba, mb) & 
!     !   +1.0d0 * t1_A(ib, bb) * t2_A(ma, nb, aa, mb) 
!        if (ia /= ma .and. aa /= na) then
!        M2_A(ib,nb,ab,bb) = M2_A(ib,nb,ab,bb) &
!        +1.0d0 * t1_A(ia, aa) * t2_A(ma, nb, ba, mb)
!        endif
!        if (ia /= ma .and. ba /= na) then
!        M2_A(ib,nb,ab,bb) = M2_A(ib,nb,ab,bb) & 
!        -1.0d0 * t1_A(ia, ba) * t2_A(ma, nb, aa, mb)
!        endif
!      enddo
!    enddo
!  enddo
!
!  !! Deltas:((nb, ib))
!      do bb = i_bb, f_bb
!        ba = bb - cc_nVa
!    do ab = i_ab, f_ab
!      aa = ab - cc_nVa
!  do jb = i_jb, f_jb
!    ja = jb - cc_nOa
!     !   M2_A(nb,jb,ab,bb) = M2_A(nb,jb,ab,bb) & 
!     !   !-1.0d0 * t1_A(jb, ab) * t2_A(ma, nb, ba, mb) & 
!     !   !+1.0d0 * t1_A(jb, bb) * t2_A(ma, nb, aa, mb) & 
!     !   !+1.0d0 * t1_A(ja, aa) * t2_A(ma, nb, ba, mb) & 
!     !   !-1.0d0 * t1_A(ja, ba) * t2_A(ma, nb, aa, mb)
!     !   +1.0d0 * t1_A(jb, ab) * t2_A(ma, nb, ba, mb) & 
!     !   -1.0d0 * t1_A(jb, bb) * t2_A(ma, nb, aa, mb)
!        if (ja /= ma .and. aa /= na) then
!        M2_A(nb,jb,ab,bb) = M2_A(nb,jb,ab,bb) & 
!        -1.0d0 * t1_A(ja, aa) * t2_A(ma, nb, ba, mb)
!        endif
!        if (ja /= ma .and. ba /= na) then
!        M2_A(nb,jb,ab,bb) = M2_A(nb,jb,ab,bb) & 
!        +1.0d0 * t1_A(ja, ba) * t2_A(ma, nb, aa, mb)
!        endif
!      enddo
!    enddo
!  enddo

end

