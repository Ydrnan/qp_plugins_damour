subroutine dm_antisym(a,lda,n,info)
        implicit none
        
        include 'constants.h'

        !======================================================
        ! Take a matrix a and give back an antisymmetric matrix
        !======================================================
        
        !===========
        ! Variables
        !===========
        
        ! in
        integer, intent(in) :: lda,n
        ! n    : integer, number of rows of the matrix a
        ! lda  : integer, specifies the leading dimension of a, must be at least max(1,n)


        !in-out 
        double precision, intent(inout) :: a(lda,n)
        ! a    : n by n double precision matrix which will be antisymmetrize
        ! in   : matrix
        ! out  : antisymmetric matrix

        ! out
        integer :: info 
        ! info : integer
                ! if info = 0, the execution is successful
                ! if info = k, the k-th parameter has an illegal value
                ! if info = -k, the algorithm failed 
        
        ! internal
        integer :: i,j
        logical :: d
        ! i,j  : integer, first and second indices to turn on the matrix elements
        ! d    : logical, display or not the matrix a after the process 

        ! external
        logical :: disnan
        ! disnan : check if an element is NaN

        !===============
        ! Pre conditions
        !===============
        
        info=0
        
        ! Size of matrix a must be at least 1 by 1
        if (n<1) then
                info = 3
                print*, 'dm_antisym : invalid parameter 3'
                print*, 'n < 1'
                return
        endif

        ! Leading dimension of A must be >= n
        if (lda < n) then
                info = 23
                print*, 'dm_antisym : invalid parameter 2 or 3'
                print*, 'lda < n'
                return
        endif
        
        ! Matrix elements must by non-NaN
        do j=1,n
                do i=1,n
                        if (disnan(a(i,j))) then 
                                info = 1
                                print*, 'dm_antisym : invalid parameter 1'
                                print*, 'NaN element in a matrix'
                                return
                        endif
                enddo
        enddo

        !============
        ! Calculation
        !============ 
        do j=1,n
                do i=1,j-1
                        a(i,j) = -a(j,i)
                enddo
                a(i,i) = 0.d0
        enddo

        !================
        ! Post conditions
        !================

        !========
        ! Display
        !========

        d=.true.
        
        if (d) then
                print*,'Matrix after antisymmetrization a :'
                do i=1,n
                        print*, a(i,:)
                enddo
        endif

end subroutine
