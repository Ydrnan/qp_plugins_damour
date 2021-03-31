program init_nb_iteration

  implicit none

  !================================================
  ! Little prog to initialize the iteration number
  ! for the trust region method and for the
  ! Umrigar method
  !================================================

  ! Start at 0
  open(unit=10,file='nb_iteration.dat')
    write(10,*) 0
  close(10)

end program
