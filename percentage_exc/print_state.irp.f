subroutine print_state(istate)

  implicit none

  BEGIN_DOC
  ! Print the first det of each state
  END_DOC
  
  ! in
  integer, intent(in) :: istate
  
  ! internal
  integer                       :: degree, h1,p1,h2,p2, exc(0:2,2)
  double precision              :: phase
  integer, allocatable          :: det(:)
  double precision, allocatable :: coef(:)
  integer                       :: i
  
  if (N_det == 1) then
    return
  endif

  allocate(det(N_det), coef(N_det))

  ! Key for dsort
  do i = 1, N_det
    det(i) = i
  enddo

  ! -abs to sort in the right order
  do i = 1, N_det
    coef(i) = -dabs(psi_coef(i,istate))
  enddo
 
  ! Sort the coef
  call dsort(coef, det, N_det)

  ! I don't like minus sign
  do i = 1, N_det
    coef(i)  = dabs(coef(i))
  enddo

  print*,''
  write(*,'(A8,I3)') '# State ', istate

  i = 1
  !print*,'psi',psi_coef(:,istate)
  do while (dabs(coef(i)) > thresh_c2) 
    print*,''
    write(*,'(A6,I10,A8,I4,A11,1pE14.6)') 'Det n. ', i, ', state n.', istate, ', |coef| = ', coef(i)

    ! Exc / ref
    call get_excitation_degree(psi_det(1,1,1),psi_det(1,1,det(i)),degree,N_int)
    write(*,'(A20,I3)') 'Excitation degree = ', degree

    ! Print exc
    if (degree <= 4 .and. degree > 0) then
      call get_excitation_spin(psi_det_alpha(1,1),psi_det_alpha(1,det(i)),exc,degree,phase,N_int) 
      call decode_exc_spin(exc,h1,p1,h2,p2)
      if (degree == 1 ) then
        write(*,'(A4,I4,A3,I4)') 'Exc:', h1, ' ->', p1
      elseif (degree == 2) then
        write(*,'(A4,I4,A3,I4,A1,I4,A3,I4)') 'Exc:', h1, ' ->', p1, ',', h2, ' ->', p2
      endif

      call get_excitation_spin(psi_det_beta(1,1),psi_det_beta(1,det(i)),exc,degree,phase,N_int)
      call decode_exc_spin(exc,h1,p1,h2,p2)
      if (degree == 1) then
        write(*,'(A4,I4,A3,I4)') 'Exc:', h1, ' ->', p1
      elseif (degree == 2) then
        write(*,'(A4,I4,A3,I4,A1,I4,A3,I4)') 'Exc:', h1, ' ->', p1, ',', h2, ' ->', p2
      endif
    endif

    ! Print det
    call print_det(psi_det(1,1,det(i)),N_int)

    i = i+1
    if (i > N_det) then
      exit 
    endif

  enddo

  print*,''

  deallocate(det, coef)

end

