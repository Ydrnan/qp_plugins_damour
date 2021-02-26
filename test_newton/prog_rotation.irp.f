program prog_rotation
      implicit none

      integer :: n,i,info
      double precision, allocatable :: A(:,:),R(:,:),C(:,:),D(:,:),B(:,:)
      double precision :: norm
      n=2

      allocate(A(n,n), B(n,n), C(n,n), D(n,n), R(n,n))


      !call random_number(A)
      A=3.14d0
      !print*,'A', A(:,:)
      call dm_antisym(A,size(A,1),n,info)      
!      if (info /= 0) then
!         call ABORT
!      endif
!      !Test 1, rotation
       
      norm = norm2(A)
      print*, 'norm :', norm


      call dm_rotation(A,size(A,1),R,size(R,1),n,info)
      call dgemm('T','N',n,n,n,1d0,R,size(R,1),R,size(R,1),0d0,B,size(B,1))
      C=0d0
      do i=1,n
        C(i,i) = 1d0
      enddo
      B=B-C
      !print*,B(:,:)
     
      do i=1,n
        write(*,'(100(F10.5))') R(i,:)
      enddo

      norm = norm2(R)
      print*, 'norm :', norm


!      if (info /= 0) then
!         call ABORT
!      endif
!      
!      !Test 2, vérification allocation matrice 
!      C=R
!      call dm_antisym(C,size(C,1),n,info)
!      if (info /= 0) then
!         call ABORT
!      endif
!      call dm_rotation(C,size(C,1),D,size(D,1),n,info)
!      if (info /= 0) then
!         call ABORT
!      endif
!
!      call dm_antisym(C,size(C,1),n,info)
!      if (info /= 0) then
!         call ABORT
!      endif
!      call dm_rotation(C,size(C,1),D,size(D,1),n,info)
!      if (info /= 0) then
!         call ABORT
!      endif

      
      deallocate(A,B,C,D,R)
      
end program
