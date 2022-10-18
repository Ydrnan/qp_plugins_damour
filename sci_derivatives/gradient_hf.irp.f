subroutine nucl_gradient_hf(dE)

  implicit none

  BEGIN_DOC
  ! Nuclear derivative of the Hatree-Fock energy
  END_DOC

  ! out
  double precision, intent(out) :: dE(nucl_num,3)

  ! internal
  double precision, allocatable :: Sa(:), Ua(:), ha(:)
  double precision, allocatable :: Ja(:,:), Ka(:,:), dE_nucl(:,:)
  integer :: i, j 
  integer :: nO, nucl, dir

  if (elec_alpha_num /= elec_beta_num) then
    print*,'Error, elec_alpha_num /= elec_beta_num, abort'
    call abort
  endif

  ! number of occupied MOs
  nO = elec_alpha_num

  ! Allocation
  allocate(Sa(nO), Ua(nO), ha(nO))
  allocate(Ja(nO,nO), Ka(nO,nO))
  allocate(dE_nucl(nucl_num,3))

  !!! Chemist notation for the bi-electronic integrals !!!

  print*,''
  print*,'The program is left as an exercise for the programmer'
  print*,''
  call abort

  ! Init
  dE = 0d0

  !!! Electronic part

  ! For each nucleus and direction (1=x, 2=y, 3=z)
  do nucl = 1, nucl_num
    do dir = 1, 3

      ! S^a
      !call 
      do i = 1, nO
        Sa(i) = 0d0 ! Sa(i,i)
      enddo
      
      ! U^a in HF
      do i = 1, nO
        Ua(i) = - 0.5d0 * Sa(i)
      enddo

      ! h_ij^a
      !call 
      do i = 1, nO
        ha(i) = 0d0 ! ha(i,i)
      enddo

      ! (ij[kl)^a
      !call

      ! (ii|jj)^a
      do i = 1, nO
        do j = 1, nO
          Ja(i,j) = 0d0 !(i,i,j,j)
        enddo
      enddo

      ! (ij|ij)^a
      do i = 1, nO
        do j = 1, nO
          Ka(i,j) = 0d0 !(i,j,i,j)
        enddo
      enddo

      ! Gradient as a 2D array of size (nucl_num x 3)
      ! dE_elec = 2 \sum_i h_ii^a + \sum_ij 2 (ii|jj)^a - (ij|ij)^a + 2 \sum_i S_ii^a e_i
      do i = 1, nO
        dE(nucl,dir) = dE(nucl,dir) + 2d0 * (ha(i) + Sa(i) * fock_matrix_diag_mo(i))
        do j = 1, nO
          dE(nucl,dir) = dE(nucl,dir) + 2d0 * Ja(i,j) - Ka(i,j)
        enddo
      enddo

    enddo
  enddo

  !!! Nuclear part
  call dE_nucl_part(dE_nucl)

  do nucl = 1, nucl_num
    do dir = 1, 3
      ! dE_nucl = \sum_{A, A /= B}  (X_A - X_B) (Z_A Z_B)/R_AB^3
      dE(nucl,dir) = dE(nucl,dir) + dE_nucl(nucl,dir)
    enddo
  enddo

  deallocate(Sa,Ua,ha,Ja,Ka)
  deallocate(dE_nucl)

end
