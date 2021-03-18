subroutine trust_region(n,method,H,v_grad,m_Hm1g)
  implicit none

  integer, intent(in) :: n 
  integer, intent(in) :: method ! pour la verif
  double precision, intent(in) :: H(n,n), v_grad(n)   
 
  double precision, intent(out) :: m_Hm1g(mo_num,mo_num)

  double precision, allocatable :: p(:),W(:,:),W_t(:,:)
  double precision, allocatable :: Hm1(:,:), gHm1(:)
  double precision :: accu, lambda, trust_radius, norm_p, norm_g, delta, rho
  double precision, allocatable :: e_val(:),work(:,:)
  integer :: nb_iter
  integer :: info,i,j,k,lwork

  ! Fonction fortran
  double precision :: ddot, norm2

  ! Fonctions 
  double precision :: fN

  lwork=3*n-1

  allocate(p(n),W(n,n),W_t(n,n))
  allocate(Hm1(n,n),gHm1(n))
  allocate(e_val(n),work(lwork,n))

  ! la diagonalisation remplace la matrice par les vecteurs propes 
  W=H
  
  ! Diagonalization
  call dsyev('V','U',n,W,size(W,1),e_val,work,lwork,info)

  if (info /= 0) then
      print*, 'Error diagonalization : trust_region'
      call ABORT
  endif

  print *, 'vp Hess:'
  write(*,'(100(F10.5))')  real(e_val(:))

  lambda =0d0

  print*,'||p||^2'
  norm_p = fN(n,e_val,W,v_grad,0d0)
  print*, norm_p

  norm_g = (norm2(v_grad))**2
  print*,'norm grad'
  print*, norm_g
 
 
  open(unit=10,file='nb_iteration.dat')
  read(10,*) nb_iter
  close(10)

  open(unit=10,file='nb_iteration.dat')
  write(10,*) nb_iter + 1
  close(10)  

  ! Qualité du modèle
  if (nb_iter == 0) then    
 
    open(unit=10,file='prev_energy.dat')
    write(10,*) 0d0
    close(10)
    
    open(unit=10,file='e_model.dat')
    write(10,*) 0d0
    close(10) 
  endif
  
  call dn_rho_model(rho)

  ! Modifs :
  ! ajouter le rejet de pas
  ! le rejet de pas doit entrainer le non changement de prev_energy
  ! Sinon l'étape suivante sera 'forcement' bonne
  ! dn_e_model doit prendre rho en paramètre 
  ! en cas d'annulation du pas il faut faire un pas qui corresponde a
  ! l'annulation du pas précédent et une division par 4 de la  
  ! region de confiance


  ! trust radius
  if (nb_iter ==0) then
    
    trust_radius = MIN(norm_g,norm_p)   

    open(unit=10,file='trust_radius.dat')
    write(10,*) trust_radius ! Delta
    close(10)

  else

    open(unit=10,file='trust_radius.dat')
    read(10,*) trust_radius ! Delta
    close(10)
    
  endif
   
  ! Passage de Delta^2 a Delta
  delta = dsqrt(trust_radius)
  print*, 'Delta', delta

  ! Changement de delta en fonction de rho
  if (rho >= 0.75) then 
    delta = 2d0 * delta
  elseif ( 0.25 >= rho ) then
    delta = delta
  else 
    delta = 0.25d0 * delta
  endif

  ! Remplacement du trust radius pour le porchaine iteration
  open(unit=10,file='trust_radius.dat')
  write(10,*) delta**2 ! Delta
  close(10)

  ! Methode de newton pour trouver le lambda pour avoir ||p|| = Delta 
  if (trust_radius < norm_p ) then 
    call trust_newton(n,e_val,W,v_grad,trust_radius,lambda)
  else 
    lambda = 0d0
  endif

  ! Initialisation
  p = 0d0 
  
  ! Calcul de p 
  do i = 1, n
    if (e_val(i) > 1d-4) then
    accu = 0d0
    accu = ddot(n,W(:,i),1,v_grad,1)
    p = p - accu * W(:,i) / (e_val(i) + lambda)
    endif
  enddo

  ! pour avoir la meme chose que gHm1
  p = -p 

  ! pour la prochaine iteration de trust region
  call dn_e_model(n,v_grad,H,p)

  ! Vector with n element -> mo_num by mo_num matrix
  do i=1,mo_num
    do j=1,mo_num
      if (i>j) then
        call in_vec_to_mat(i,j,k)
        m_Hm1g(i,j) = p(k)
      else
        m_Hm1g(i,j)=0d0
      endif
    enddo
  enddo

  ! Antisymmetrization of the previous matrix
  do i=1,mo_num
    do j=1,mo_num
      if (i<j) then
        m_Hm1g(i,j) = - m_Hm1g(j,i)
      endif
    enddo
  enddo

! Debug
  print*,'p'
  write(*,'(100(F10.5))') p(:)
  
  ! Verification
  call dm_inversion(method,n,H,Hm1)

  print*,''
  call dgemv('T',n,n,1d0,Hm1,size(Hm1,1),v_grad,1,0d0,gHm1,1)

  print*,'vector g^T.Hm1 :'
  write(*,'(100(F10.5))') gHm1(:)
  !print*, gHm1

  ! Calcul de la différence / erreur
  p = p - gHm1
  print*,'diff'
  do i = 1, n
    if (ABS(p(i)) > 1e-12) then
      print*,i,p(i)
    endif
  enddo

end subroutine
