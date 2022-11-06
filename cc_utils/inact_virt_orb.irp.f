BEGIN_PROVIDER [integer, n_inact_virt_no_core_orb ]
  implicit none
  BEGIN_DOC
  !  Number of inactive and virtual MOs
  END_DOC
  n_inact_virt_no_core_orb = (n_inact_orb + n_virt_orb)
END_PROVIDER 

BEGIN_PROVIDER [integer, dim_list_inact_virt_no_core_orb]
   implicit none
   BEGIN_DOC
   ! dimensions for the allocation of list_inact_virt_no_core
   ! it is at least 1
   END_DOC
   dim_list_inact_virt_no_core_orb = max(n_inact_virt_no_core_orb,1)
END_PROVIDER

BEGIN_PROVIDER [integer, list_inact_virt_no_core, (dim_list_inact_virt_no_core_orb)]
 implicit none
 integer :: i,j
 j = 0
 do i = 1, mo_num
  if(  trim(mo_class(i))=="Inactive" &
  .or. trim(mo_class(i))=="Virtual" )then
   j += 1
   list_inact_virt_no_core(j) = i
  endif
 enddo

END_PROVIDER
