program init_cyrus
  implicit none
   
  integer :: n,i
  n= mo_num*(mo_num-1)/2

  open(unit=10,file='prev2_Hm1g.dat')
  do i=1,n
    write(10,*) 0.0d0
  enddo
  close(10)

  open(unit=11,file='prev_Hm1g.dat')
  do i=1,n
    write(11,*) 0.0d0
  enddo
  close(11)
end program
