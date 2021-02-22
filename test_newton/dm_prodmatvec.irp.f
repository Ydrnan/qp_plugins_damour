subroutine dm_prodvecmat(v,n,A,LDA,B,LDB,top,info)
        implicit none

        !===========
        ! Variables
        !===========

        ! int
        integer :: n, LDA, LDB
        double precision, intent(in)  :: v(n)
        ! n   : integer, number of element of the vector that contain the
        !       diagonal elements of diagonal matrix/ eigenvalues
        ! LDA : integer, specifies the leading dimension of A, must be at least max(1,n) 
        ! v   : double precision vector with the diagonal elements of 
        !       diagonal matrix/ eigenvalue

        ! inout       
        double precision, intent(inout) :: A(LDA,n)
        ! A   : n by n double precision matrix which will be multiply by v
        !       and store the result
        
        !out
        integer, intent(out) :: info
        double precision, intent(out) :: B(LDB,n)
        ! info : integer :
                ! if info = 0, the execution is successful
                ! if info = k, the k-th parameter has an illegal value
                ! if info = -k, the algorithm failed
 
        ! internal
        !double precision, allocatable :: B(:,:)
        integer :: i,j,top
        logical :: d
        ! B   : n by n double precision matrix, use to compute the
        !       the product with transpose matrix without compute 
        !       the transpose matrix of A
        ! i,j : integer, first and second indices to turn on the matrix elements
        ! top : integer, type operation 
        !       - 0 : diagonal matrix . matrix
        !       - 1 : diagonal matrix . matrix^t
        !       - 2 : matrix . diagonal matrix
        !       - 3 : matrix . diagonal matrix^t
        ! d   : logical, if true display the rotation matrix
        
        ! external
        logical :: disnan
        ! disnan : check if an element is NaN 

        !===============
        ! Preconditions
        !===============

        ! Initialization
        info = 0 

        ! Size of matrix A must be at least 1 by 1
        if (n < 1) then
                info = 2
                print*, 'dm_prodmatvec : invalid argument 2'
                print*, 'n < 1'
                return
        endif

        ! Leading dimension of A must be >= n
        if (LDA < n) then 
                info = 24
                print*, 'dm_prodmatvec : invalid argument 2 or 4'
                print*, 'LDA < n'
                return
        endif

        ! Leading dimension of B must be >= n
        if (LDB < n) then
                info = 26
                print*, 'dm_prodmatvec : invalid argument 2 or 6'
                print*, 'LDB < n'
                return
        endif

        ! Vector elements of v must be non-Nan
        do i=1,n
                if (disnan(v(i))) then
                        info = 1
                        print*, 'dm_prodmatvec : invalid parameter 1'
                        print*, 'NaN element in v vector'
                        return
                endif
        enddo

        ! Matrix elements of A must by non-NaN
        do j=1,n
                do i=1,n
                        if (disnan(A(i,j))) then
                                info = 3
                                print*, 'dm_prodmatvec : invalid parameter 3'
                                print*, 'NaN element in A matrix'
                                return
                        endif
                enddo
        enddo
        
        !==============
        ! Calculations
        !==============

        if (top==0) then ! diagonal matrix . matrix
                do j=1,n
                         do i=1,n
                                B(i,j)=v(i)*A(i,j)
                         enddo
                enddo
                !do i=1,n
                !        call dscal(n, v(i), A(:,i), 1)
                !enddo

        else if (top==1) then ! diagonal matrix . matrix^t
                do j=1,n
                         do i=1,n
                                B(i,j)=v(i)*A(j,i)
                         enddo
                enddo
        
        elseif (top==2) then ! matrix . diagonal matrix
                do j=1,n
                         do i=1,n
                                B(i,j)=A(i,j)*v(j)
                         enddo
                enddo

        elseif (top==3) then ! matrix . diagonal matrix^t
                do j=1,n
                         do i=1,n
                                B(i,j)=A(i,j)*v(i)
                         enddo
                enddo
        else
                info = 5
                print*, 'dm_prodmatvec : Unknown argument 7'
        endif


        !=================
        ! Post conditions
        !=================

        !=========
        ! Display
        !=========
        
        if (d) then 
                print*, 'A'
                do i=1,n
                        print*, A(i,:)
                enddo
        endif
end subroutine
