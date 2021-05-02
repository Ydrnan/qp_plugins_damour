program extrapolation
implicit none
integer :: nb_lines, i, k, nb_points, index_first_PT2
double precision, allocatable :: E(:), E_corr(:), PT2(:), rPT2(:)
character(len=10) :: text_PT2, text_rPT2
character(len=100) :: ezfio_dir
character(len=200) :: path_cipsi
character(len=120) :: path_hf
double precision :: hf_e

!write(*,*) 'ezfio directory name ?'
!read(*,*) ezfio_dir
!print*, ezfio_dir

!path_hf = ezfio_dir//'hartree_fock/energy'

! ! on utilise une variable logique
 logical::existe
!
! ! on teste l'existence de mon_fichier
! inquire( file=path_hf, exist=existe)
!
! ! la réponse est dans la variable existe
! If    ( existe ) Then
!    ! existe = .true.
!    write(*,*)"le répertoire existe"
! ElseIf( .NOT. existe) Then
!    ! existe = .false.
!    write(*,*)"le répertoire n'existe pas"
!    STOP 
! Endif

write(*,*) 'fci.out name ?'
read(*,*) path_cipsi

!path_cipsi = adjustl(adjustr(ezfio_dir)//'/'//path_cipsi)

!print*, path_cipsi
call execute_command_line('echo '//path_cipsi)
 ! on teste l'existence de mon_fichier
 inquire( file=path_cipsi, exist=existe)

 ! la réponse est dans la variable existe
 If    ( existe ) Then
    ! existe = .true.
    write(*,*)"le fichier existe"
 ElseIf( .NOT. existe) Then
    ! existe = .false.
    write(*,*)"le fichier n'existe pas"
    STOP
 Endif

call execute_command_line('grep "# E  " '//path_cipsi//' > E.dat')
!call execute_command_line("grep '# E  ' CN_6_31g.fci.out > E.dat")
call execute_command_line("awk '{print $3}' E.dat > newE.dat")

call execute_command_line('grep "# PT2  " '//path_cipsi//' > PT2.dat')
!call execute_command_line("grep '# PT2  ' CN_6_31g.fci.out > PT2.dat")
call execute_command_line("awk '{print $3}' PT2.dat > newPT2.dat")

call execute_command_line('grep "# rPT2  " '//path_cipsi//' > rPT2.dat')
!call execute_command_line("grep '# rPT2  ' CN_6_31g.fci.out > rPT2.dat")
call execute_command_line("awk '{print $3}' rPT2.dat > newrPT2.dat")

call execute_command_line('grep "Summary at N_det = " '//path_cipsi//' > Ndet.dat')
!call execute_command_line("grep 'Summary at N_det = ' CN_6_31g.fci.out > Ndet.dat")
call execute_command_line("awk '{print $5}' Ndet.dat > newNdet.dat")

!call execute_command_line("sed 's/#//' E.dat > tmp1.dat")
call execute_command_line("wc -l E.dat > nb_lines.dat")
!call execute_command_line("sed 's/tmp1.dat//' nb_line.dat > example2.dat")
call execute_command_line("awk '{print $1}' nb_lines.dat > newnb_lines.dat")

! Number of lines
open(unit=10,file='newnb_lines.dat')
read(10,*) nb_lines
close(10)
print*,nb_lines

! allocation
allocate(E(nb_lines))
allocate(PT2(nb_lines))
allocate(rPT2(nb_lines))
allocate(E_corr(nb_lines))

! read the ci energies
open(unit=10,file='newE.dat')
do i = 1, nb_lines
  read(10,*) E(i)
enddo
close(10)
print*,'E'
print*,E(:)

open(unit=10,file='newPT2.dat')
do i = 1, nb_lines
  read(10,*) PT2(i)
enddo
close(10)
print*,'PT2'
print*,PT2(:)

open(unit=10,file='newrPT2.dat')
do i = 1, nb_lines
  read(10,*) rPT2(i)
enddo
close(10)
print*,'rPT2'
print*,rPT2(:)

!hf_e = hf_energy
hf_e = -88d0
!path_hf = ezfio_dir//'hartree_fock/energy'
!open(unit=10,file=path_hf)
!read(10,*) hf_e
!close(10)

! compute and write the correlation energies (mH)
do i = 1, nb_lines
  E_corr(i) = (E(i) - hf_e) * 1d3
enddo
print*,'e_corr'
print*,E_corr(:)

open(unit=10,file='e_correlation.dat')
do i = 1, nb_lines
  write(10,*) E_corr(i)
enddo
close(10)

call execute_command_line("paste newNdet.dat newE.dat > tmp2.dat")
call execute_command_line("paste tmp2.dat newPT2.dat > tmp3.dat")
call execute_command_line("paste tmp3.dat newrPT2.dat > tmp4.dat")
call execute_command_line("paste tmp4.dat e_correlation.dat > tmp5.dat")
call execute_command_line("cat newnb_lines.dat > data_correlation.dat")

open(unit=10,file='data_correlation.dat')
write(10,*) nb_lines
write(10,*) "N_det      E           PT2          rPT2           E_corr"
close(10)

call execute_command_line("cat tmp5.dat >> data_correlation.dat")

call execute_command_line("rm Ndet.dat E.dat PT2.dat rPT2.dat nb_lines.dat")
call execute_command_line("rm newNdet.dat newE.dat newPT2.dat newrPT2.dat e_correlation.dat newnb_lines.dat")
call execute_command_line("rm tmp2.dat tmp3.dat tmp4.dat tmp5.dat")

k = MAX(1, nb_lines-7)

open(unit=10, file='E_correlation.dat')
write(10,*) 'FCI energy extrapolation'
write(10,*) ''
write(10,*) 'With PT2 :'
write(10,*) 'Nb points  E_corr     R^2     error    error^2'
write(10,*) '            (mH)              (mH)     (mH^2)'
close(10)

index_first_PT2 = MAX(1,nb_lines-6)

write(text_PT2,'(F10.4)') PT2(index_first_PT2)
write(text_rPT2,'(F10.4)') rPT2(index_first_PT2)

open(unit=10, file='E_correlation.gnu')
write(10,*) '#!/bin/gnuplot'
write(10,*) 'set xrange ['//text_PT2//':0]'
write(10,*) 'plot "data_correlation.dat" u 3:5 title "E corr = f(PT2)"'
close(10)

do nb_points = 3, nb_lines - k + 1
call linear_reg(nb_lines - k + 1 - nb_points + 1, nb_lines, PT2, E_corr)
enddo

open(unit=10, file='E_correlation.dat', position='append')
write(10,*) ''
write(10,*) 'With rPT2 :'
write(10,*) 'Nb points  E_corr     R^2     error    error^2'
write(10,*) '            (mH)              (mH)     (mH^2)'
close(10)

do nb_points = 3, nb_lines - k + 1
call linear_reg(nb_lines - k + 1 - nb_points + 1, nb_lines, rPT2, E_corr)
enddo

end program

subroutine linear_reg(k,nb_lines,x,y)
implicit none

  ! Linear regression to fin a and b such as :
  ! y = a*x + b

  ! a = \sum_i (x_i - \bar{x})(y_i - \bar{y}) / \sum_i (x_i - \bar{x})^2
  ! b = \bar_y - a \bar{x}

  ! \bar{x} = 1/ n \sum_i x_i 
  ! \bar{y} = 1/ n \sum_i y_i 

  ! R^2 = 1 - SSR / SST
  
  ! SSR = \sum_i (y_i - \hat{y})^2
  ! SST = \sum_i (y_i - \bar{y})^2

  ! \hat{y}_i = a*x_i + b

  integer, intent(in) :: k, nb_lines
  double precision, intent(in) :: x(nb_lines), y(nb_lines)

  double precision :: a,b
  double precision :: x_bar, y_bar, num, denom
  double precision :: SSR, SST, R2, abs_error, abs_error2
  double precision, allocatable :: y_hat(:)
  integer :: i, n
  character(len=10) :: text_a, text_b
  character(len=2)  :: text_points
  character(len=10) :: text_title1
  character(len=3) :: text_title2
  character(len=10) :: text_title3
  character(len=14) :: text

  ! Number of points for the regression
  n = nb_lines - k + 1
  print*,n
  ! Allocation
  allocate(y_hat(nb_lines))

  ! x_bar
  x_bar = 0d0
  do i = k, nb_lines
    x_bar = x_bar + x(i)
  enddo
  x_bar = x_bar / DBLE(n)

  ! y_bar
  y_bar = 0d0
  do i = k, nb_lines
    y_bar = y_bar + y(i)
  enddo
  y_bar = y_bar / DBLE(n)

  ! a = num/denom
  num = 0d0
  do i = k, nb_lines
    num = num + (x(i)-x_bar)*(y(i)-y_bar)
  enddo

  denom = 0d0
  do i = k, nb_lines
    denom = denom + (x(i)-x_bar)**2
  enddo

  a = num/denom

  b = y_bar - a * x_bar

  ! Check
  print*,'a :',a
  print*,'b :',b 

  do i = k, nb_lines
    print*,'Reality :', y(i), 'Model :', a*x(i) + b
  enddo

  ! y_hat
  y_hat = 0d0
  do i = k, nb_lines
    y_hat(i) = a*x(i) + b
  enddo
  
  ! abs_error
  ! Sum of the absolute error between the reality and the model
  ! divided by the number of points
  abs_error = 0d0
  do i = k, nb_lines
    abs_error = abs_error + ABS(y(i) - y_hat(i))
  enddo

  abs_error = abs_error / DBLE(n)

  ! SSR
  ! Sum of square residuals
  SSR = 0d0
  do i = k, nb_lines
    SSR = SSR + (y(i) - y_hat(i))**2
  enddo

  ! abs_error2
  ! Sum of the error^2 between the reality and the model
  ! divided by the number of points
  abs_error2 = SSR / DBLE(n)
  
  ! SST
  ! Sum of square total
  SST = 0d0
  do i = k, nb_lines
    SST = SST + (y(i) - y_bar)**2
  enddo

  ! R^2, quality of the model
  R2 = 1d0 - SSR/SST

  ! Display
  print*,'R^2 :', R2
  print*,'Error (mH) :', abs_error
  print*,'Error^2 :', abs_error2

  open(unit=10, file='E_correlation.dat', position='append')
  write(10,'(I7, F10.2, F10.5, F10.5, F10.5)') n, b, R2, abs_error, abs_error2              
  close(10)

  write(text_a,'(F10.3)') a 
  write(text_b,'(F10.3)') b
  write(text_points,'(I2)') n 
  text_points = trim(text_points)
  text_title1 = ' title "'
  text_title2 = '"'
  write(text_title3,'(F10.2)') b
  text = adjustl(adjustr(text_title3)//' mH'//'"')
  print*, text

  open(unit=10, file='E_correlation.gnu', position='append')
  write(10,*) 'replot ', text_a, '*x', text_b, text_title1, text_points,' points : ', text
  close(10)

end subroutine
