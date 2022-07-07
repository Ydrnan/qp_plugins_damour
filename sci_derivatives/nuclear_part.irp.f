subroutine dE_nucl_part(dE_nucl)

  implicit none

  BEGIN_DOC
  ! Derivative of the nuclear repulsion energy
  END_DOC

  ! inout
  double precision, intent(inout) :: dE_nucl(nucl_num,3)

  ! internal
  integer :: nucl, dir, A, B
  double precision :: dir_AB, r_AB

  dE_nucl = 0d0

  ! dE_nuc = \sum_{A, A /= B}  (dir_A - dir_B) (Z_A Z_B)/R_AB^3
  do dir = 1, 3
    do A = 1, nucl_num
      do B = 1, nucl_num
        if (A == B) then
          cycle
        endif
        dir_AB = nucl_coord(A,dir) - nucl_coord(B,dir)
        r_AB = nucl_dist(A,B)

        dE_nucl(A,dir) = dE_nucl(A,dir) - dir_AB * (nucl_charge(A) * nucl_charge(B))/r_AB**3
      enddo
    enddo
  enddo

end
