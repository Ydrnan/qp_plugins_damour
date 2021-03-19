program init_nb_iteration
  implicit none

   open(unit=10,file='nb_iteration.dat')
     write(10,*) 0
   close(10)

end program
