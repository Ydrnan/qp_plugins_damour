subroutine trust_region(n,method,H,v_grad,m_Hm1g)

  include 'constants.h' 

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

  if (debug) then
    print *, 'vp Hess:'
    write(*,'(100(F10.5))')  real(e_val(:))
  endif

  lambda =0d0

  print*,'||p||^2 :'
  norm_p = fN(n,e_val,W,v_grad,0d0)
  print*, norm_p

  norm_g = (norm2(v_grad))**2
  print*,'norm grad :'
  print*, norm_g
 
  open(unit=10,file='nb_iteration.dat')
    read(10,*) nb_iter
  close(10)

  open(unit=10,file='nb_iteration.dat')
    if (nb_iter == -1) then
      write(10,*) nb_iter + 2
    else
      write(10,*) nb_iter + 1
    endif
  close(10)  

  ! Qualité du modèle
  call dn_rho_model(rho,nb_iter)

  ! trust radius
  if (nb_iter ==0) then
    
    trust_radius = norm_p !MIN(norm_g,norm_p)   

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
  print*, 'delta', delta

  ! Changement de delta en fonction de rho
  if (rho >= 0.75d0) then 
    delta = 2d0 * delta
  elseif (rho >= 0.5d0) then
    delta = delta
  elseif (rho >= 0.25d0) then
    delta = 0.5d0 * delta
  else
    delta = 0.25d0 * delta
  endif

  ! Remplacement du trust radius pour le porchaine iteration
  open(unit=10,file='trust_radius.dat')
    write(10,*) delta**2 ! Delta
  close(10)

  ! En donnant delta, on cherche ||p||^2 - delta^2 = 0
  ! et non ||p||^2 - delta = 0
  print*,'trust_radius',trust_radius
  print*,'delta * coef',delta

  if (rho >= 0.1d0) then ! rho minimum pour accepter les pas !0.25d0
 
    ! Methode de newton pour trouver le lambda pour avoir ||p|| = Delta 
    if (trust_radius < norm_p ) then 
      !call trust_newton(n,e_val,W,v_grad,trust_radius,lambda)
      call trust_newton(n,e_val,W,v_grad,delta,lambda)
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
    !p = 10d0 * p ! test annulation de pas

    ! Ecriture du pas
    open(unit=10,file='Hm1g.dat')
      do i = 1, n
        write(10,*) p(i)
      enddo
    close(10)
  
  else 
   
    ! Lecture du dernier pas pour l'annuler
    open(unit=10,file='Hm1g.dat')
      do i = 1, n
        read(10,*) p(i)
      enddo
    close(10)
  
    ! Annulation du pas dernier pas
    !p = -0.75d0 * p 
    p = -p 

    ! Remplacement de e_model et prev en simulant un premier pas 
    ! pour la prochaine itération
 
    open(unit=10,file='nb_iteration.dat')
      write(10,*) -1
    close(10)

  endif

  ! Energie prevue par le modele pour la prochaine iteration 
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
  if (debug) then 
 
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
 
 endif

end subroutine
