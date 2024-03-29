subroutine rotation_matrix_omp(A,LDA,R,LDR,n,info)
     
        use omp_lib     
   
        implicit none

        include 'constants.h'
        
        !======================================================
        ! Build a rotation matrix from an antisymmetric matrix
        !======================================================

        !===========
        ! Variables
        !===========

        ! in
        integer, intent(in)           :: n,LDA,LDR
        ! n   : integer, number of columns of the A matrix
        ! LDA : integer, specifies the leading dimension of A, must be at least max(1,n)
        ! LDR : integer, specifies the leading dimension of R, must be at least max(1,n)

        ! out
        double precision, intent(out) :: A(LDA,n), R(LDR,n)
        integer, intent(out)          :: info
        ! A    : n by n antisymmetric double precision matrix
        ! R    : n by n rotation double precision matrix
        ! info : integer :
                ! if info = 0, the execution is successful
                ! if info = k, the k-th parameter has an illegal value
                ! if info = -k, the algorithm failed

        ! internal
        double precision, allocatable :: B(:,:) ! initial matrix A, B=A.A
        double precision, allocatable :: work(:,:) ! work matrix for dsyev
        double precision, allocatable :: W(:,:), e_val(:)!,W_t(:,:) ! Eigenvectors and eigenvalues
        double precision, allocatable :: m_diag(:,:),cos_tau(:,:),sin_tau(:,:),tau_m1(:,:)
        double precision, allocatable :: v_cos_tau(:), v_sin_tau(:), v_taum1(:)
        double precision, allocatable :: part_1(:,:),part_1a(:,:)
        double precision, allocatable :: part_2(:,:),part_2a(:,:),part_2b(:,:),part_2c(:,:)
        double precision, allocatable :: RR_t(:,:)
        double precision, allocatable :: WW_t(:,:) ! pour test
        double precision, dimension(1:2,1:2) :: rot_2 ! test
        integer :: i,j
        integer :: info2, lwork ! for dsyev
        double precision :: t1,t2
        double precision :: norm, norm_2

        ! B       : n by n symmetric double precision matrix, B=A.A
        ! lwork   : integer, dimension of the syev work array >= max(1, 3n-1)
        ! work    : double precision matrix for syev, dimension max(1,lwork)
        ! W       : n by n double precision matrix with the eigenvectors of B diagonalization
        !
        ! e_val   : n, double precision vector with the eigenvalues of B diagonalization
        ! m_diag  : n by n double precision diagonal matrix of B diagonalization
        ! cos_tau : n by n double precision diagonal matrix with cos(tau) values
        ! sin_tau : n by n double precision diagonal matrix with sin(tau) values
        ! tau_m1  : n by n double precision diagonal matrix with 1/tau values
        ! part_1  : n by n double precision matrix, W.cos_tau.W^t
        ! part_1a : n by n double precision matrix, cos_tau.W^t
        ! part_2  : n by n double precision matrix, W.tau_m1.sin_tau.W^t.A
        ! part_2a : n by n double precision matrix, W^t.A
        ! part_2b : n by n double precision matrix, sin_tau.W^t.A
        ! part_2c : n by n double precision matrix, tau_m1.sin_tau.W^t.A
        ! RR_t    : n by n double precision matrix, R.R^t must be equal to the identity
        !           if R is a rotation matrix <=> R.R^t-1=0 <=> norm = 0
        ! norm    : integer, norm of R.R^t-1, must be equal to 0
        ! i,j     : integer, first and second indices to turn on the matrix elements

        ! intrinsic
        double precision :: dnrm2
        logical :: disnan
        ! dnrm2 : double precision function, compute the norm of a matrix
        ! disnan : check if an element is NaN

        !============
        ! Allocation
        !============

        allocate(B(n,n))
        allocate(m_diag(n,n),cos_tau(n,n),sin_tau(n,n),tau_m1(n,n))
        allocate(v_cos_tau(n),v_sin_tau(n),v_taum1(n))
        allocate(W(n,n),WW_t(n,n))
        allocate(part_1(n,n),part_1a(n,n))
        allocate(part_2(n,n),part_2a(n,n),part_2b(n,n),part_2c(n,n))
        allocate(RR_t(n,n))

        !================
        ! Pre-conditions
        !================

        if (debug) then
          print*,'Enter in rotation_matrix_omp'
        endif

        info=0

        ! Size of matrix A must be at least 1 by 1
        if (n<1) then
                info = 3
                print*, 'rotation_matrix_omp : invalid parameter 5'
                print*, 'n<1'
                return
        endif

        ! Leading dimension of A must be >= n
        if (LDA < n) then
                info = 25
                print*, 'rotation_matrix_omp : invalid parameter 2 or 5'
                print*, 'LDA < n'
                return
        endif

        ! Leading dimension of A must be >= n
        if (LDR < n) then
                info = 4
                print*, 'rotation_matrix_omp : invalid parameter 4'
                print*, 'LDR < n'
                return
        endif

        ! Matrix elements of A must by non-NaN
        call omp_set_max_active_levels(1)

        !$OMP PARALLEL                                                     &
          !$OMP PRIVATE(i,j)                                               &
          !$OMP SHARED(a,info,n)&
          !$OMP DEFAULT(NONE)

        !$OMP DO
        do j=1,n
          do i=1,n
            if (disnan(a(i,j))) then
              info=1
              print*, 'rotation_matrix_omp : invalid parameter 1'
              print*, 'NaN element in A matrix'
            endif
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
       
        if (info == 1) then
          return
        endif

        !==============
        ! Calculations
        !==============

        !==============
        ! Compute B=A.A
        
        !$OMP PARALLEL     &
          !$OMP SHARED(n,A,B) &
          !$OMP DEFAULT(NONE)
        call dgemm('N','N',n,n,n,1d0,A,size(A,1),A,size(A,1),0d0,B,size(B,1))
        !$OMP END PARALLEL
     

        ! Debug
        ! Display the B=A.A
        if (debug) then
          print*,'B=A.A'
          do i=1,n
            print*, B(i,:)
          enddo
        endif

        !============================================================
        ! Copy B in W, diagonalization will put the eigenvectors in W
        W=B

        !=====================
        ! Diagonalization of B
        ! Eigenvalues -> e_val
        ! Eigenvectors -> W
        lwork=3*n-1
        allocate(work(lwork,n),e_val(n))

        print*,'Starting diagonalization ...'

        call dsyev('V','U',n,W,size(W,1),e_val,work,lwork,info2)
           
        deallocate(work)

        if (info2==0) then
          print*, 'Diagonalization : Done'
        elseif (info<0) then
          print*, 'Diagonalization : error'
          print*, 'Illegal value of the ', info2,'-th parameter'
        else
          print*, "Diagonalization : Failed to converge"
        endif

        ! Debug
        if (debug) then
          print*, 'Eigenvalues'
          print*, e_val(:)
          print*, 'Eigenvectors'
          do i=1,n
            print*, W(i,:)
          enddo
        endif

        ! W.W^t-1 =? 0
        if (debug) then
          do j=1,n
            do i=1,n
              if (i==j) then
                WW_t(i,j)=1d0
              else
                WW_t(i,j)=0d0
              endif
            enddo
          enddo
          
          call dgemm('N','T',n,n,n,1d0,W,size(W,1),W,size(W,1),-1d0,WW_t,size(WW_t,1))
          norm = dnrm2(n*n,WW_t,1) / (dble(n)**2)
          print*, 'norm = ', norm

          print*, 'W.W^t'
          do i=1,n
                  print*,WW_t(i,:)
          enddo
        endif 

        ! ======================
        ! Diagonal matrix m_diag

        do j=1,n
          if (e_val(j) >= 0.d0) then
            e_val(j) = 0.d0
          else
             e_val(j) = -e_val(j)
          endif
        enddo

        m_diag = 0.d0
        do i=1,n
          m_diag(i,i)=e_val(i)
        enddo

        !=================================================
        ! eigenvalues = -tau^2 => tau = -sqrt(eignevalues)
        ! We need the diagonal matrix cos_tau and sin_tau, such as :
        ! if i==j then cos_tau(i,j)=cos(tau(i)) else 0d0
        ! if i==j then sin_tau(i,j)=sin(tau(i)) else 0d0

        !$OMP PARALLEL     &
          !$OMP PRIVATE(i,j) &
          !$OMP SHARED(n,cos_tau,e_val) &
          !$OMP DEFAULT(NONE)

        !$OMP DO
        do j=1,n
          do i=1,n
            if (i==j) then
              cos_tau(i,j)=dcos(dsqrt(e_val(i)))
            else
              cos_tau(i,j)=0d0
            endif
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        do i=1,n
          v_cos_tau(i)=dcos(dsqrt(e_val(i)))
        enddo

        !$OMP PARALLEL     &
          !$OMP PRIVATE(i,j) &
          !$OMP SHARED(n,sin_tau,e_val) &
          !$OMP DEFAULT(NONE)

        !$OMP DO
        do j=1,n
          do i=1,n
            if (i==j) then
              sin_tau(i,j)=dsin(dsqrt(e_val(i)))
            else
              sin_tau(i,j)=0d0
            endif
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        do i=1,n
          v_sin_tau(i)=dsin(dsqrt(e_val(i)))
        enddo

        ! Debug
        ! Display the cos_tau and sin_tau matrix
        if (debug) then
          print*, 'cos_tau'
          do i=1,n
                  print*, cos_tau(i,:)
          enddo
          print*, 'sin_tau'
          do i=1,n
                  print*, sin_tau(i,:)
          enddo
        endif

        !=======
        ! tau^-1

        !$OMP PARALLEL     &
          !$OMP PRIVATE(i,j) &
          !$OMP SHARED(n,tau_m1,e_val) &
          !$OMP DEFAULT(NONE)

        !$OMP DO
        do j=1,n
          do i=1,n
            if ((i==j).and.(e_val(i) > 0.d0)) then
              tau_m1(i,j)=1d0/(dsqrt(e_val(i)))
            else
              tau_m1(i,j)=0d0
            endif
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        do i=1,n
          if (e_val(i) > 0.d0) then
            v_taum1(i)=1d0/(dsqrt(e_val(i)))
          else
            v_taum1(i)=0d0
          endif
        enddo

        ! Debug
        ! Display tau^-1
        if (debug) then
          print*, 'tau^-1'
          do i=1,n
            print*,tau_m1(i,:)
          enddo
        endif

        !============================================================
        ! We now have to compute W.cos_tau.W^t+W.tau_m1.sin_tau.W^t.A
        ! First part_1=dgemm(W, dgemm(cos_tau, W^t))
        ! part_1a = dgemm(cos_tau, W^t)
        ! part_1 = dgemm(W, part_1a)
        ! And part_2=dgemm(W, dgemm(tau_m1, dgemm(sin_tau, dgemm(W^t, A))))
        ! part_2a = dgemm(W^t, A)
        ! part_2b = dgemm(sin_tau, part_2a)
        ! part_2c = dgemm(tau_m1, part_2b)
        ! part_2 = dgemm(W, part_2c)
        ! Rotation matrix R = part_1+part_2

        !call cpu_time(t1)
        !$OMP PARALLEL     &
          !$OMP SHARED(n,cos_tau,sin_tau,tau_m1,W,part_1a,part_1, &
          !$OMP part_2a,part_2b,part_2c,part_2,A) &
          !$OMP DEFAULT(NONE) 
        call dgemm('N','T',n,n,n,1d0,cos_tau,size(cos_tau,1),W,size(W,1),0d0,part_1a,size(part_1a,1))
        call dgemm('N','N',n,n,n,1d0,W,size(W,1),part_1a,size(part_1a,1),0d0,part_1,size(part_1,1))

        call dgemm('T','N',n,n,n,1d0,W,size(W,1),A,size(A,1),0d0,part_2a,size(part_2a,1))
        call dgemm('N','N',n,n,n,1d0,sin_tau,size(sin_tau,1),part_2a,size(part_2a,1),0d0,part_2b,size(part_2b,1))
        call dgemm('N','N',n,n,n,1d0,tau_m1,size(tau_m1,1),part_2b,size(part_2b,1),0d0,part_2c,size(part_2c,1))
        call dgemm('N','N',n,n,n,1d0,W,size(W,1),part_2c,size(part_2c,1),0d0,part_2,size(part_2,1))
        !$OMP END PARALLEL

        !=========================
        ! n by n rotation matrix R
        R = part_1 + part_2

        !=============
        ! Matrix check
        ! R.R^t and R^t.R must be equal to identity matrix
        ! On ne fait la vérification que sur R.R^t pour le moment

        !$OMP PARALLEL     &
          !$OMP PRIVATE(i,j) &
          !$OMP SHARED(n,RR_t,R) &
          !$OMP DEFAULT(NONE)

        !$OMP DO
        do j=1,n
          do i=1,n
            if (i==j) then
              RR_t(i,j)=1d0
            else
              RR_t(i,j)=0d0
            endif
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
  
        call dgemm('N','T',n,n,n,1d0,R,size(R,1),R,size(R,1),-1d0,RR_t,size(RR_t,1))

        norm = dnrm2(n*n,RR_t,1) / (dble(n)**2)
        print*, 'Rotation matrix check, norm R.R^T = ', norm

        ! Debug
        if (debug) then
          print*, 'RR_t'
          do i=1,n
            print*, RR_t(i,:)
          enddo
        endif

        !=================
        ! Post-conditions
        !=================

        ! Matrix elements of R must by non-NaN

        !$OMP PARALLEL     &
          !$OMP PRIVATE(i,j) &
          !$OMP SHARED(n,LDR,R,info) &
          !$OMP DEFAULT(NONE)

        !$OMP DO
        do j=1,n
          do i=1,LDR
            if (disnan(R(i,j))) then
              info=666
              print*, 'NaN in rotation matrix'
              call ABORT
            endif
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        !=========
        ! Display
        !=========

        if (debug) then
          print*,'Rotation matrix :'
          do i=1,mo_num
            write(*,'(100(F10.5))') R(i,:)
          enddo
        endif

        !==============
        ! Deallocation
        !==============

        deallocate(B)
        deallocate(m_diag,cos_tau,sin_tau,tau_m1)
        deallocate(v_cos_tau,v_sin_tau,v_taum1)
        deallocate(W,WW_t)
        deallocate(part_1,part_1a)
        deallocate(part_2,part_2a,part_2b,part_2c)
        deallocate(RR_t)

        if (debug) then
          print*,'Leave rotation_matrix_omp'
        endif

end subroutine dm_rotation

