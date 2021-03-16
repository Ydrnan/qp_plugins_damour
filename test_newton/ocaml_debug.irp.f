subroutine ocaml_debug()

  include 'constants.h'
  
  implicit none
 
  integer :: i,j

  if (ocaml) then 
    open(unit=10,file='../../../../../../App_y/miniconda3/Work_yann/mo_num.dat')
      write(10,*) mo_num
    close(10)
  
    open(unit=11,file='../../../../../../App_y/miniconda3/Work_yann/mo_coef.dat')
    do i = 1, mo_num
      do j = 1, mo_num
        write(11,*) i,j, mo_coef(i,j)
      enddo
    enddo
    close(11)
  endif
end subroutine
