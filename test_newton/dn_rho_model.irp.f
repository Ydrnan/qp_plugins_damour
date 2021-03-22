subroutine dn_rho_model(rho,nb_iter)
  implicit none
  !calcul le rapport : (prev_energy - energy) / (prev_energy - e_modele) 
  
  double precision, intent(out) :: rho

  integer, intent(in) :: nb_iter
  
  double precision :: energy, prev_energy, e_model

  !open(unit=10,file='energy.dat')
  !read(10,*) energy
  !close(10)
 
  energy = ci_energy(1) 
  print*, 'ci_energy', energy

  if ( nb_iter >=1) then 
  
    open(unit=11,file='prev_energy.dat')
    read(11,*) prev_energy
    close(11)
  
    open(unit=12,file='e_model.dat')
    read(12,*) e_model
    close(12)
    
    rho = (prev_energy - energy) / (prev_energy - e_model)
    
    print*, 'prev_energy :', prev_energy
    print*, 'e_model :', e_model
    print*, 'energy :', energy
    print*, 'prev_energy - energy', prev_energy - energy
    print*, 'e_model - energy', prev_energy - e_model
    print*,'Rho :', rho
  
  elseif (nb_iter == -1) then 
    rho = 0.5d0

    open(unit=11,file='prev_energy.dat')
      read(11,*) prev_energy
    close(11)

    print*, 'prev_energy :', prev_energy
  else
    rho = 0.5d0
  endif

  if (rho >= 0.1d0 .and. nb_iter >= 0) then

    open(unit=11,file='prev_energy.dat')
      write(11,*) energy
      print*, 'energie -> prev_energy :', energy
    close(11)
  endif
   
end subroutine
