subroutine dn_e_model(n,v_grad,H,Hm1g)
  implicit none

  integer, intent(in) :: n
  double precision, intent(in) :: v_grad(n),H(n,n),Hm1g(n)

  double precision :: e_model, prev_energy
  double precision :: ddot
  double precision :: part_1, part_2
  double precision, allocatable :: part_2a(:)

  integer :: i,j

  allocate(part_2a(n))

  !!!!!! Il faut automatiser l'écriture de l'énergie fci pour que cela marche !!!!!!
  open(unit=10,file='prev_energy.dat')
    read(10,*) prev_energy
  close(10)

  part_1 = ddot(n,v_grad,1,Hm1g,1)
 
 !print*,'v_grad'
  !print*, v_grad(:)
  !print*,'Hm1g'
  !print*, Hm1g(:)

  print*,'part_1 : ', part_1
  
  call dgemv('N',n,n,1d0,H,size(H,1),Hm1g,1,0d0,part_2a,1)
  
  part_2 = 0.5d0 * ddot(n,Hm1g,1,part_2a,1)
  
  print*,'part_2 : ', part_2 

  ! Verif, pourquoi part_1 et 2 sont positifs ??
  ! A revoir
  e_model = prev_energy - part_1 - part_2

  print*, 'e_model : ', e_model

  open(unit=10,file='e_model.dat')
    write(10,*) e_model
  close(10)

 deallocate(part_2a)

end subroutine 
