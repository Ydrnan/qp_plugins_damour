subroutine dm_vec_to_mat(v,k,B,n,info)
        implicit none

        !=============================================================
        ! Build a n by n antisymmetric matrix from a n*(n-1)/2 vector 
        !=============================================================

        !===========
        ! Variables
        !===========

        ! in
        integer, intent(in)           :: k,n
        double precision, intent(in)  :: v(k)
        ! k : integer, number of elements in the vector v
        ! n : size of the matrix create from v

        ! out 
        double precision, intent(out) :: B(n,n)
        integer, intent(out)          :: info
        ! B    : n by n double precision antisymmetric matrix from v
        ! info : integer
                ! if info = 0, the execution is successful
                ! if info = k, the k-th parameter has an illegal value

        ! internal
        integer :: i,j,l
        logical :: d
        ! i,j : integer, first and second indices to turn on the matrix elements
        ! l   : integer, to scan the vector v
        ! d   : logical, if true display the matrix

        ! external
        logical :: disnan
        ! disnan : check if an element is NaN

        !================
        ! Pre conditions
        !================
        
        ! Size of vector a must be at least 1 
        if (k<1) then
                info = 2
                print*, 'dm_vec_mat : invalid parameter 2'
                print*, 'k < 1'
                return
        endif

        ! Size of vector must be k=n*(n-1)/2
        if (k/=(n*(n-1)/2)) then
                info = 24
                print*, 'dm_vec_mat : invalid parameter 2 or 4'
                print*, 'k/=(n*(n-1)/2)'
                return
        endif

        ! Vector elements must by non-NaN
        do i=1,k
                if (disnan(v(k))) then
                        info = 1
                        print*, 'dm_vec_mat : invalid parameter 1'
                        print*, 'NaN element in v vector'
                        return
                endif
        enddo
       
        !=============
        !Calculations
        !=============

        ! Vector -> lower diagonal matrix
        l=1
        do j=1,n-1
                do i=j+1,n
                        B(i,j)=v(l)
                        l=l+1
                enddo
        enddo

        ! Lower diagonal matrix -> antisymmetric matrix 
        do j=1,n
                do i=1,j-1
                        B(i,j) = -B(j,i)
                enddo
                B(j,j) = 0.d0
        enddo

        !=================
        ! Post conditions
        !=================

        !=========
        ! Display
        !=========
        
        d=.true.

        !if (d) then
        !        print*,'B'
        !        do i=1,n
        !                print*,B(i,:)
        !        enddo
        !endif

end subroutine
