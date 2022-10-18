
! ---

double precision function overlap_1d_G_Gxx(A_center, B_center, alpha, beta, power_A, power_B0, dim)

  BEGIN_DOC
  !
  ! \int_{-infty}^{+infty} Gx(A) \partial_x^2[ Gx(B) ] dx
  !
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,          intent(in) :: dim, power_A, power_B0
  double precision, intent(in) :: A_center, B_center, alpha, beta

  integer                      :: i, iorder_p, power_B
  double precision             :: P_new(0:max_dim), P_center, fact_p, p
  double precision             :: integ_tmp

  double precision             :: F_integral


  overlap_1d_G_Gxx = 0.d0
  power_B = power_B0

  ! ---

  call give_explicit_poly_and_gaussian_x(P_new, P_center, p, fact_p, iorder_p, alpha, beta, power_A, power_B, A_center, B_center, dim)

  if(fact_p .lt. 1.d-15) then
    overlap_1d_G_Gxx = 0.d0
    return
  endif

  integ_tmp = 0.d0
  do i = 0, iorder_p
    integ_tmp = integ_tmp + P_new(i) * F_integral(i, p)
  enddo
  overlap_1d_G_Gxx = overlap_1d_G_Gxx & 
                   + ( -2.d0 * beta * (2.d0*dble(power_B) + 1.d0) ) * integ_tmp * fact_p

  ! ---

  if(power_B .gt. -1) then

    power_B = power_B - 2
    call give_explicit_poly_and_gaussian_x(P_new, P_center, p, fact_p, iorder_p, alpha, beta, power_A, power_B, A_center, B_center, dim)
    power_B = power_B + 2

    integ_tmp = 0.d0
    do i = 0, iorder_p
      integ_tmp = integ_tmp + P_new(i) * F_integral(i, p)
    enddo
    overlap_1d_G_Gxx = overlap_1d_G_Gxx & 
                     + ( dble(power_B) * (dble(power_B) - 1.d0) ) * integ_tmp * fact_p

  endif

  ! ---

  power_B = power_B + 2
  call give_explicit_poly_and_gaussian_x(P_new, P_center, p, fact_p, iorder_p, alpha, beta, power_A, power_B, A_center, B_center, dim)
  power_B = power_B - 2

  integ_tmp = 0.d0
  do i = 0, iorder_p
    integ_tmp = integ_tmp + P_new(i) * F_integral(i, p)
  enddo
  overlap_1d_G_Gxx = overlap_1d_G_Gxx & 
                   + ( 4.d0 * beta * beta ) * integ_tmp * fact_p

  ! ---

end function overlap_1d_G_Gxx

! ---

subroutine overlap_3d_G_Gxx(A_center, B_center, alpha, beta, power_A, power_B, dim, Ix, Iy, Iz)

  BEGIN_DOC
  !
  !   Ix = \int_{-infty}^{+infty} Gx(A) \partial_x^2[ Gx(B) ] dx x 
  !        \int_{-infty}^{+infty} Gy(A) Gy(B) dy                 x
  !        \int_{-infty}^{+infty} Gz(A) Gz(B) dz
  !
  !   Iy = \int_{-infty}^{+infty} Gx(A) Gx(B) dx                 x 
  !        \int_{-infty}^{+infty} Gy(A) \partial_y^2[ Gy(B) ] dy x
  !        \int_{-infty}^{+infty} Gz(A) Gz(B) dz
  !
  !   Iz = \int_{-infty}^{+infty} Gx(A) Gx(B) dx                 x 
  !        \int_{-infty}^{+infty} Gy(A) Gy(B) dy                 x
  !        \int_{-infty}^{+infty} Gz(A) \partial_z^2[ Gz(B) ] dz
  !
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,          intent(in)  :: dim, power_A(3), power_B(3) 
  double precision, intent(in)  :: A_center(3), B_center(3), alpha, beta
  double precision, intent(out) :: Ix, Iy, Iz

  integer                      :: i, nmax, iorder_p(3), iorder_p_x
  integer                      :: power_A_x, power_B_x
  double precision             :: A_center_x, B_center_x, P_center_x
  double precision             :: P_new(0:max_dim,3), P_center(3), P_new_x(0:max_dim), fact_p, p
  double precision             :: F_integral_tab(0:max_dim)
  double precision             :: integ_tmp
  double precision             :: Ix0, Iy0, Iz0
  double precision             :: Ix1, Iy1, Iz1
  double precision             :: Ixm, Iym, Izm
  double precision             :: Ixp, Iyp, Izp

  double precision             :: F_integral


  Ix = 0.d0
  Iy = 0.d0
  Iz = 0.d0

  call give_explicit_poly_and_gaussian(P_new, P_center, p, fact_p, iorder_p, alpha, beta, power_A, power_B, A_center, B_center, dim)
  if(fact_p .lt. 1d-15) then
     return
  endif

  ! p = alpha + beta for all the 1d-integrals
  F_integral_tab = 0.d0
  nmax = maxval(iorder_p) + 2
  do i = 0, nmax
    F_integral_tab(i) = F_integral(i, p)
  enddo

  ! ---

  ! x term 

  call gaussian_product_x(alpha, A_center(1), beta, B_center(1), fact_p, p, P_center(1))
  integ_tmp = 0.d0
  do i = 0, iorder_p(1)
    integ_tmp = integ_tmp + P_new(i,1) * F_integral_tab(i)
  enddo
  Ix0 = integ_tmp * fact_p
  Ix1 = Ix0 * ( -2.d0 * beta * (2.d0*dble(power_B(1)) + 1.d0) ) 

  ! y term 

  call gaussian_product_x(alpha, A_center(2), beta, B_center(2), fact_p, p, P_center(2))
  integ_tmp = 0.d0
  do i = 0, iorder_p(2)
    integ_tmp = integ_tmp + P_new(i,2) * F_integral_tab(i)
  enddo
  Iy0 = integ_tmp * fact_p
  Iy1 = Iy0 * ( -2.d0 * beta * (2.d0*dble(power_B(2)) + 1.d0) ) 

  ! z term 

  call gaussian_product_x(alpha, A_center(3), beta, B_center(3), fact_p, p, P_center(3))
  integ_tmp = 0.d0
  do i = 0, iorder_p(3)
    integ_tmp = integ_tmp + P_new(i,3) * F_integral_tab(i)
  enddo
  Iz0 = integ_tmp * fact_p
  Iz1 = Iz0 * ( -2.d0 * beta * (2.d0*dble(power_B(3)) + 1.d0) ) 

  ! ---

  power_A_x  = power_A (1)
  power_B_x  = power_B (1)
  A_center_x = A_center(1)
  B_center_x = B_center(1)

  if(power_B_x .gt. -1) then

    power_B_x = power_B_x - 2
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 2

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo
    Ixm = ( dble(power_B_x) * (dble(power_B_x) - 1.d0) ) * integ_tmp * fact_p

  else
 
    Ixm = 0.d0

  endif

  power_B_x = power_B_x + 2
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 2

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo
  Ixp = ( 4.d0 * beta * beta ) * integ_tmp * fact_p

  ! ---

  power_A_x  = power_A (2)
  power_B_x  = power_B (2)
  A_center_x = A_center(2)
  B_center_x = B_center(2)

  if(power_B_x .gt. -1) then

    power_B_x = power_B_x - 2
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 2

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo
    Iym = ( dble(power_B_x) * (dble(power_B_x) - 1.d0) ) * integ_tmp * fact_p

  else
 
    Iym = 0.d0

  endif

  power_B_x = power_B_x + 2
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 2

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo
  Iyp = ( 4.d0 * beta * beta ) * integ_tmp * fact_p

  ! ---

  power_A_x  = power_A (3)
  power_B_x  = power_B (3)
  A_center_x = A_center(3)
  B_center_x = B_center(3)

  if(power_B_x .gt. -1) then

    power_B_x = power_B_x - 2
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 2

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo
    Izm = ( dble(power_B_x) * (dble(power_B_x) - 1.d0) ) * integ_tmp * fact_p

  else
 
    Izm = 0.d0

  endif

  power_B_x = power_B_x + 2
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 2

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo
  Izp = ( 4.d0 * beta * beta ) * integ_tmp * fact_p

  ! ---

  Ix = (Ixm + Ix1 + Ixp) * Iy0 * Iz0
  Iy = (Iym + Iy1 + Iyp) * Ix0 * Iz0
  Iz = (Izm + Iz1 + Izp) * Iy0 * Ix0

end subroutine overlap_3d_G_Gxx

! ---

subroutine overlap_3d_G_GXbxx(A_center, B_center, alpha, beta, power_A, power_B, dim, Ix, Iy, Iz)

  BEGIN_DOC
  !
  !   Ix = \int_{-infty}^{+infty} Gx(A) \partial_{XB} \partial_x^2[ Gx(B) ] dx x 
  !        \int_{-infty}^{+infty} Gy(A) Gy(B) dy                               x
  !        \int_{-infty}^{+infty} Gz(A) Gz(B) dz
  !
  !   Iy = \int_{-infty}^{+infty} Gx(A) \partial_{XB} Gx(B) dx   x 
  !        \int_{-infty}^{+infty} Gy(A) \partial_y^2[ Gy(B) ] dy x
  !        \int_{-infty}^{+infty} Gz(A) Gz(B) dz
  !
  !   Iz = \int_{-infty}^{+infty} Gx(A) \partial_{XB} Gx(B) dx x 
  !        \int_{-infty}^{+infty} Gy(A) Gy(B) dy               x
  !        \int_{-infty}^{+infty} Gz(A) \partial_z^2[ Gz(B) ] dz
  !
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,          intent(in)  :: dim, power_A(3), power_B(3) 
  double precision, intent(in)  :: A_center(3), B_center(3), alpha, beta
  double precision, intent(out) :: Ix, Iy, Iz

  integer                      :: i, nmax, iorder_p(3), iorder_p_x
  integer                      :: power_A_x, power_B_x
  double precision             :: A_center_x, B_center_x, P_center_x
  double precision             :: P_new(0:max_dim,3), P_center(3), P_new_x(0:max_dim), fact_p, p
  double precision             :: F_integral_tab(0:max_dim)
  double precision             :: integ_tmp
  double precision             :: Ix0, Iy0, Iz0
  double precision             :: Ix1, Iy1, Iz1
  double precision             :: Ixm3, Iym3, Izm3
  double precision             :: Ixm2, Iym2, Izm2
  double precision             :: Ixm1, Iym1, Izm1
  double precision             :: Ixp1, Iyp1, Izp1
  double precision             :: Ixp2, Iyp2, Izp2
  double precision             :: Ixp3, Iyp3, Izp3

  double precision             :: F_integral


  Ix = 0.d0
  Iy = 0.d0
  Iz = 0.d0

  call give_explicit_poly_and_gaussian(P_new, P_center, p, fact_p, iorder_p, alpha, beta, power_A, power_B, A_center, B_center, dim)
  if(fact_p .lt. 1d-15) then
     return
  endif

  ! p = alpha + beta for all the 1d-integrals
  F_integral_tab = 0.d0
  nmax = maxval(iorder_p) + 3
  do i = 0, nmax
    F_integral_tab(i) = F_integral(i, p)
  enddo

  ! ---

  call gaussian_product_x(alpha, A_center(1), beta, B_center(1), fact_p, p, P_center(1))

  integ_tmp = 0.d0
  do i = 0, iorder_p(1)
    integ_tmp = integ_tmp + P_new(i,1) * F_integral_tab(i)
  enddo

  Ix0 = integ_tmp * fact_p

  ! ---

  call gaussian_product_x(alpha, A_center(2), beta, B_center(2), fact_p, p, P_center(2))

  integ_tmp = 0.d0
  do i = 0, iorder_p(2)
    integ_tmp = integ_tmp + P_new(i,2) * F_integral_tab(i)
  enddo

  Iy0 = integ_tmp * fact_p
  Iy1 = Iy0 * ( -2.d0 * beta * (2.d0*dble(power_B(2)) + 1.d0) ) 

  ! ---

  call gaussian_product_x(alpha, A_center(3), beta, B_center(3), fact_p, p, P_center(3))

  integ_tmp = 0.d0
  do i = 0, iorder_p(3)
    integ_tmp = integ_tmp + P_new(i,3) * F_integral_tab(i)
  enddo

  Iz0 = integ_tmp * fact_p
  Iz1 = Iz0 * ( -2.d0 * beta * (2.d0*dble(power_B(3)) + 1.d0) ) 

  ! ---

  ! ---------------------------------------------------------------------------------------------------------
  ! d_x^2

  power_A_x  = power_A (1)
  power_B_x  = power_B (1)
  A_center_x = A_center(1)
  B_center_x = B_center(1)

  ! ---

  if(power_B_x .gt. -2) then

    power_B_x = power_B_x - 3
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 3

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo

    Ixm3 = ( -dble(power_B_x) * (dble(power_B_x) - 1.d0) * (dble(power_B_x) - 2.d0) ) * integ_tmp * fact_p

  else
 
    Ixm3 = 0.d0

  endif

  ! ---

  if(power_B_x .gt. 0) then

    power_B_x = power_B_x - 1
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 1

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo

    Ixm1 = ( 6.d0 * beta * dble(power_B_x) * dble(power_B_x) ) * integ_tmp * fact_p
    Iym1 = ( -dble(power_B_x) ) * integ_tmp * fact_p
    Izm1 = Iym1

  else
 
    Ixm1 = 0.d0
    Iym1 = 0.d0
    Izm1 = 0.d0

  endif

  ! ---

  power_B_x = power_B_x + 1
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 1

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo

  Ixp1 = ( -12.d0 * beta * beta * (dble(power_B_x) + 1.d0) ) * integ_tmp * fact_p
  Iyp1 = 2.d0 * beta * integ_tmp * fact_p
  Izp1 = Iyp1

  ! ---

  power_B_x = power_B_x + 3
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 3

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo

  Ixp3 = ( 8.d0 * beta * beta * beta ) * integ_tmp * fact_p

  ! ---------------------------------------------------------------------------------------------------------

  ! ---

  ! ---------------------------------------------------------------------------------------------------------
  ! d_y^2

  power_A_x  = power_A (2)
  power_B_x  = power_B (2)
  A_center_x = A_center(2)
  B_center_x = B_center(2)

  ! ---

  if(power_B_x .gt. -1) then

    power_B_x = power_B_x - 2
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 2

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo

    Iym2 = ( dble(power_B_x) * (dble(power_B_x) - 1.d0) ) * integ_tmp * fact_p

  else
 
    Iym2 = 0.d0

  endif

  ! ---

  power_B_x = power_B_x + 2
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 2

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo

  Iyp2 = 4.d0 * beta * beta * integ_tmp * fact_p

  ! ---------------------------------------------------------------------------------------------------------

  ! ---

  ! ---------------------------------------------------------------------------------------------------------
  ! d_z^2

  power_A_x  = power_A (3)
  power_B_x  = power_B (3)
  A_center_x = A_center(3)
  B_center_x = B_center(3)

  if(power_B_x .gt. -1) then

    power_B_x = power_B_x - 2
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 2

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo

    Izm2 = ( dble(power_B_x) * (dble(power_B_x) - 1.d0) ) * integ_tmp * fact_p

  else
 
    Izm2 = 0.d0

  endif

  ! ---

  power_B_x = power_B_x + 2
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 2

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo

  Izp2 = 4.d0 * beta * beta * integ_tmp * fact_p

  ! ---------------------------------------------------------------------------------------------------------

  Ix = (Ixm3 + Ixm1 + Ixp1 + Ixp3) * Iy0                 * Iz0
  Iy = (Iym1 + Iyp1)               * (Iym2 + Iy1 + Iyp2) * Iz0
  Iz = (Izm1 + Izp1)               * Iy0                 * (Izm2 + Iz1 + Izp2)

end subroutine overlap_3d_G_GXbxx

! ---

subroutine overlap_3d_G_GYbxx(A_center, B_center, alpha, beta, power_A, power_B, dim, Ix, Iy, Iz)

  BEGIN_DOC
  !
  !   Ix = \int_{-infty}^{+infty} Gx(A) \partial_x^2[ Gx(B) ] dx x 
  !        \int_{-infty}^{+infty} Gy(A) \partial_{YB} Gy(B) dy   x
  !        \int_{-infty}^{+infty} Gz(A) Gz(B) dz
  !
  !   Iy = \int_{-infty}^{+infty} Gx(A) Gx(B) dx                               x 
  !        \int_{-infty}^{+infty} Gy(A) \partial_{YB} \partial_y^2[ Gy(B) ] dy x
  !        \int_{-infty}^{+infty} Gz(A) Gz(B) dz
  !
  !   Iz = \int_{-infty}^{+infty} Gx(A) Gx(B) dx               x 
  !        \int_{-infty}^{+infty} Gy(A) \partial_{YB} Gy(B) dy x
  !        \int_{-infty}^{+infty} Gz(A) \partial_z^2[ Gz(B) ] dz
  !
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,          intent(in)  :: dim, power_A(3), power_B(3) 
  double precision, intent(in)  :: A_center(3), B_center(3), alpha, beta
  double precision, intent(out) :: Ix, Iy, Iz

  integer                      :: i, nmax, iorder_p(3), iorder_p_x
  integer                      :: power_A_x, power_B_x
  double precision             :: A_center_x, B_center_x, P_center_x
  double precision             :: P_new(0:max_dim,3), P_center(3), P_new_x(0:max_dim), fact_p, p
  double precision             :: F_integral_tab(0:max_dim)
  double precision             :: integ_tmp
  double precision             :: Ix0, Iy0, Iz0
  double precision             :: Ix1, Iy1, Iz1
  double precision             :: Ixm3, Iym3, Izm3
  double precision             :: Ixm2, Iym2, Izm2
  double precision             :: Ixm1, Iym1, Izm1
  double precision             :: Ixp1, Iyp1, Izp1
  double precision             :: Ixp2, Iyp2, Izp2
  double precision             :: Ixp3, Iyp3, Izp3

  double precision             :: F_integral


  Ix = 0.d0
  Iy = 0.d0
  Iz = 0.d0

  call give_explicit_poly_and_gaussian(P_new, P_center, p, fact_p, iorder_p, alpha, beta, power_A, power_B, A_center, B_center, dim)
  if(fact_p .lt. 1d-15) then
     return
  endif

  ! p = alpha + beta for all the 1d-integrals
  F_integral_tab = 0.d0
  nmax = maxval(iorder_p) + 3
  do i = 0, nmax
    F_integral_tab(i) = F_integral(i, p)
  enddo

  ! ---

  call gaussian_product_x(alpha, A_center(1), beta, B_center(1), fact_p, p, P_center(1))

  integ_tmp = 0.d0
  do i = 0, iorder_p(1)
    integ_tmp = integ_tmp + P_new(i,1) * F_integral_tab(i)
  enddo

  Ix0 = integ_tmp * fact_p
  Ix1 = Ix0 * ( -2.d0 * beta * (2.d0*dble(power_B(1)) + 1.d0) ) 

  ! ---

  call gaussian_product_x(alpha, A_center(2), beta, B_center(2), fact_p, p, P_center(2))

  integ_tmp = 0.d0
  do i = 0, iorder_p(2)
    integ_tmp = integ_tmp + P_new(i,2) * F_integral_tab(i)
  enddo

  Iy0 = integ_tmp * fact_p

  ! ---

  call gaussian_product_x(alpha, A_center(3), beta, B_center(3), fact_p, p, P_center(3))

  integ_tmp = 0.d0
  do i = 0, iorder_p(3)
    integ_tmp = integ_tmp + P_new(i,3) * F_integral_tab(i)
  enddo

  Iz0 = integ_tmp * fact_p
  Iz1 = Iz0 * ( -2.d0 * beta * (2.d0*dble(power_B(3)) + 1.d0) ) 

  ! ---

  ! ---------------------------------------------------------------------------------------------------------
  ! d_y^2

  power_A_x  = power_A (2)
  power_B_x  = power_B (2)
  A_center_x = A_center(2)
  B_center_x = B_center(2)

  ! ---

  if(power_B_x .gt. -2) then

    power_B_x = power_B_x - 3
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 3

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo

    Iym3 = ( -dble(power_B_x) * (dble(power_B_x) - 1.d0) * (dble(power_B_x) - 2.d0) ) * integ_tmp * fact_p

  else
 
    Iym3 = 0.d0

  endif

  ! ---

  if(power_B_x .gt. 0) then

    power_B_x = power_B_x - 1
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 1

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo

    Ixm1 = ( -dble(power_B_x) ) * integ_tmp * fact_p
    Iym1 = ( 6.d0 * beta * dble(power_B_x) * dble(power_B_x) ) * integ_tmp * fact_p
    Izm1 = Ixm1

  else
 
    Ixm1 = 0.d0
    Iym1 = 0.d0
    Izm1 = 0.d0

  endif

  ! ---

  power_B_x = power_B_x + 1
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 1

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo

  Ixp1 = 2.d0 * beta * integ_tmp * fact_p
  Iyp1 = ( -12.d0 * beta * beta * (dble(power_B_x) + 1.d0) ) * integ_tmp * fact_p
  Izp1 = Ixp1

  ! ---

  power_B_x = power_B_x + 3
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 3

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo

  Iyp3 = ( 8.d0 * beta * beta * beta ) * integ_tmp * fact_p

  ! ---------------------------------------------------------------------------------------------------------

  ! ---

  ! ---------------------------------------------------------------------------------------------------------
  ! d_x^2

  power_A_x  = power_A (1)
  power_B_x  = power_B (1)
  A_center_x = A_center(1)
  B_center_x = B_center(1)

  if(power_B_x .gt. -1) then

    power_B_x = power_B_x - 2
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 2

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo

    Ixm2 = ( dble(power_B_x) * (dble(power_B_x) - 1.d0) ) * integ_tmp * fact_p

  else
 
    Ixm2 = 0.d0

  endif

  ! ---

  power_B_x = power_B_x + 2
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 2

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo

  Ixp2 = 4.d0 * beta * beta * integ_tmp * fact_p

  ! ---------------------------------------------------------------------------------------------------------

  ! ---

  ! ---------------------------------------------------------------------------------------------------------
  ! d_z^2

  power_A_x  = power_A (3)
  power_B_x  = power_B (3)
  A_center_x = A_center(3)
  B_center_x = B_center(3)

  if(power_B_x .gt. -1) then

    power_B_x = power_B_x - 2
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 2

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo

    Izm2 = ( dble(power_B_x) * (dble(power_B_x) - 1.d0) ) * integ_tmp * fact_p

  else
 
    Izm2 = 0.d0

  endif

  ! ---

  power_B_x = power_B_x + 2
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 2

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo

  Izp2 = 4.d0 * beta * beta * integ_tmp * fact_p

  ! ---------------------------------------------------------------------------------------------------------

  Ix = (Ixm2 + Ix1 + Ixp2) * (Ixm1 + Ixp1)               * Iz0
  Iy = Ix0                 * (Iym3 + Iym1 + Iyp1 + Iyp3) * Iz0
  Iz = Ix0                 * (Izm1 + Izp1)               * (Izm2 + Iz1 + Izp2)

end subroutine overlap_3d_G_GYbxx

! ---

subroutine overlap_3d_G_GZbxx(A_center, B_center, alpha, beta, power_A, power_B, dim, Ix, Iy, Iz)

  BEGIN_DOC
  !
  !   Ix = \int_{-infty}^{+infty} Gx(A) \partial_x^2[ Gx(B) ] dx x 
  !        \int_{-infty}^{+infty} Gy(A) \partial_{YB} Gy(B) dy   x
  !        \int_{-infty}^{+infty} Gz(A) Gz(B) dz
  !
  !   Iy = \int_{-infty}^{+infty} Gx(A) Gx(B) dx                               x 
  !        \int_{-infty}^{+infty} Gy(A) \partial_{YB} \partial_y^2[ Gy(B) ] dy x
  !        \int_{-infty}^{+infty} Gz(A) Gz(B) dz
  !
  !   Iz = \int_{-infty}^{+infty} Gx(A) Gx(B) dx               x 
  !        \int_{-infty}^{+infty} Gy(A) \partial_{YB} Gy(B) dy x
  !        \int_{-infty}^{+infty} Gz(A) \partial_z^2[ Gz(B) ] dz
  !
  END_DOC

  implicit none
  include 'constants.include.F'

  integer,          intent(in)  :: dim, power_A(3), power_B(3) 
  double precision, intent(in)  :: A_center(3), B_center(3), alpha, beta
  double precision, intent(out) :: Ix, Iy, Iz

  integer                      :: i, nmax, iorder_p(3), iorder_p_x
  integer                      :: power_A_x, power_B_x
  double precision             :: A_center_x, B_center_x, P_center_x
  double precision             :: P_new(0:max_dim,3), P_center(3), P_new_x(0:max_dim), fact_p, p
  double precision             :: F_integral_tab(0:max_dim)
  double precision             :: integ_tmp
  double precision             :: Ix0, Iy0, Iz0
  double precision             :: Ix1, Iy1, Iz1
  double precision             :: Ixm3, Iym3, Izm3
  double precision             :: Ixm2, Iym2, Izm2
  double precision             :: Ixm1, Iym1, Izm1
  double precision             :: Ixp1, Iyp1, Izp1
  double precision             :: Ixp2, Iyp2, Izp2
  double precision             :: Ixp3, Iyp3, Izp3

  double precision             :: F_integral


  Ix = 0.d0
  Iy = 0.d0
  Iz = 0.d0

  call give_explicit_poly_and_gaussian(P_new, P_center, p, fact_p, iorder_p, alpha, beta, power_A, power_B, A_center, B_center, dim)
  if(fact_p .lt. 1d-15) then
     return
  endif

  ! p = alpha + beta for all the 1d-integrals
  F_integral_tab = 0.d0
  nmax = maxval(iorder_p) + 3
  do i = 0, nmax
    F_integral_tab(i) = F_integral(i, p)
  enddo

  ! ---

  call gaussian_product_x(alpha, A_center(1), beta, B_center(1), fact_p, p, P_center(1))

  integ_tmp = 0.d0
  do i = 0, iorder_p(1)
    integ_tmp = integ_tmp + P_new(i,1) * F_integral_tab(i)
  enddo

  Ix0 = integ_tmp * fact_p
  Ix1 = Ix0 * ( -2.d0 * beta * (2.d0*dble(power_B(1)) + 1.d0) ) 

  ! ---

  call gaussian_product_x(alpha, A_center(2), beta, B_center(2), fact_p, p, P_center(2))

  integ_tmp = 0.d0
  do i = 0, iorder_p(2)
    integ_tmp = integ_tmp + P_new(i,2) * F_integral_tab(i)
  enddo

  Iy0 = integ_tmp * fact_p
  Iy1 = Iy0 * ( -2.d0 * beta * (2.d0*dble(power_B(2)) + 1.d0) ) 

  ! ---

  call gaussian_product_x(alpha, A_center(3), beta, B_center(3), fact_p, p, P_center(3))

  integ_tmp = 0.d0
  do i = 0, iorder_p(3)
    integ_tmp = integ_tmp + P_new(i,3) * F_integral_tab(i)
  enddo

  Iz0 = integ_tmp * fact_p

  ! ---

  ! ---------------------------------------------------------------------------------------------------------
  ! d_z^2

  power_A_x  = power_A (3)
  power_B_x  = power_B (3)
  A_center_x = A_center(3)
  B_center_x = B_center(3)

  ! ---

  if(power_B_x .gt. -2) then

    power_B_x = power_B_x - 3
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 3

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo

    Izm3 = ( -dble(power_B_x) * (dble(power_B_x) - 1.d0) * (dble(power_B_x) - 2.d0) ) * integ_tmp * fact_p

  else
 
    Izm3 = 0.d0

  endif

  ! ---

  if(power_B_x .gt. 0) then

    power_B_x = power_B_x - 1
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 1

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo

    Ixm1 = ( -dble(power_B_x) ) * integ_tmp * fact_p
    Iym1 = Ixm1
    Izm1 = ( 6.d0 * beta * dble(power_B_x) * dble(power_B_x) ) * integ_tmp * fact_p

  else
 
    Ixm1 = 0.d0
    Iym1 = 0.d0
    Izm1 = 0.d0

  endif

  ! ---

  power_B_x = power_B_x + 1
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 1

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo

  Ixp1 = 2.d0 * beta * integ_tmp * fact_p
  Iyp1 = Ixp1
  Izp1 = ( -12.d0 * beta * beta * (dble(power_B_x) + 1.d0) ) * integ_tmp * fact_p

  ! ---

  power_B_x = power_B_x + 3
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 3

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo

  Izp3 = ( 8.d0 * beta * beta * beta ) * integ_tmp * fact_p

  ! ---------------------------------------------------------------------------------------------------------

  ! ---

  ! ---------------------------------------------------------------------------------------------------------
  ! d_x^2

  power_A_x  = power_A (1)
  power_B_x  = power_B (1)
  A_center_x = A_center(1)
  B_center_x = B_center(1)

  if(power_B_x .gt. -1) then

    power_B_x = power_B_x - 2
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim )
    power_B_x = power_B_x + 2

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo

    Ixm2 = ( dble(power_B_x) * (dble(power_B_x) - 1.d0) ) * integ_tmp * fact_p

  else
 
    Ixm2 = 0.d0

  endif

  ! ---

  power_B_x = power_B_x + 2
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 2

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo

  Ixp2 = 4.d0 * beta * beta * integ_tmp * fact_p


  ! ---------------------------------------------------------------------------------------------------------

  ! ---

  ! ---------------------------------------------------------------------------------------------------------
  ! d_y^2

  power_A_x  = power_A (2)
  power_B_x  = power_B (2)
  A_center_x = A_center(2)
  B_center_x = B_center(2)

  if(power_B_x .gt. -1) then

    power_B_x = power_B_x - 2
    call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                          , power_A_x, power_B_x, A_center_x, B_center_x, dim)
    power_B_x = power_B_x + 2

    integ_tmp = 0.d0
    do i = 0, iorder_p_x
      integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
    enddo

    Iym2 = ( dble(power_B_x) * (dble(power_B_x) - 1.d0) ) * integ_tmp * fact_p

  else
 
    Iym2 = 0.d0

  endif

  ! ---

  power_B_x = power_B_x + 2
  call give_explicit_poly_and_gaussian_x( P_new_x, P_center_x, p, fact_p, iorder_p_x, alpha, beta &
                                        , power_A_x, power_B_x, A_center_x, B_center_x, dim)
  power_B_x = power_B_x - 2

  integ_tmp = 0.d0
  do i = 0, iorder_p_x
    integ_tmp = integ_tmp + P_new_x(i) * F_integral_tab(i)
  enddo

  Iyp2 = 4.d0 * beta * beta * integ_tmp * fact_p

  ! ---------------------------------------------------------------------------------------------------------

  Ix = (Ixm2 + Ix1 + Ixp2) * Iy0                 * (Ixm1 + Ixp1)
  Iy = Ix0                 * (Iym2 + Iy1 + Iyp2) * (Iym1 + Iyp1)
  Iz = Ix0                 * Iy0                 * (Izm3 + Izm1 + Izp1 + Izp3)

end subroutine overlap_3d_G_GZbxx

! ---


