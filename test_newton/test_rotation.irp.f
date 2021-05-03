program test_rotation
      implicit none

      integer :: n,i,info
      double precision, allocatable :: A(:,:),R(:,:),B(:,:)
      double precision :: norm, dnrm2
      n=2

      allocate(A(n,n), B(n,n), R(n,n))


      A=3.14d0
      call dm_antisym(A,size(A,1),n,info)             
      call rotation_matrix(A,size(A,1),R,size(R,1),n,info)

      ! R^T.R = identity
      call dgemm('T','N',n,n,n,1d0,R,size(R,1),R,size(R,1),0d0,B,size(B,1))
      do i = 1, n
        B(i,i) = B(i,i) - 1d0
      enddo
      norm = dnrm2(n,B,1)
      print*,'Norm R^T.R:', norm    
     
       print*,'R :' 
      do i=1,n
        write(*,'(100(F10.5))') R(i,:)
      enddo

      deallocate(A,B,R)
      
end program
