program toto
      implicit none

      integer :: n,info
      double precision, allocatable :: A(:,:),R(:,:),C(:,:),D(:,:),B(:,:)

      n=3

      allocate(A(n,n), B(n,n), C(n,n), D(n,n), R(n,n))


      call random_number(A)
      print*,'A', A(:,:)
      call dm_antisym(A,size(A,1),n,info)      
      if (info /= 0) then
         call ABORT
      endif
      !Test 1, rotation
      call dm_rotation(A,size(A,1),R,size(R,1),n,info)
      if (info /= 0) then
         call ABORT
      endif
      
      !Test 2, vérification allocation matrice 
      C=R
      call dm_antisym(C,size(C,1),n,info)
      if (info /= 0) then
         call ABORT
      endif
      call dm_rotation(C,size(C,1),D,size(D,1),n,info)
      if (info /= 0) then
         call ABORT
      endif

      call dm_antisym(C,size(C,1),n,info)
      if (info /= 0) then
         call ABORT
      endif
      call dm_rotation(C,size(C,1),D,size(D,1),n,info)
      if (info /= 0) then
         call ABORT
      endif

      
      deallocate(A,B,C,D,R)
      
end program
