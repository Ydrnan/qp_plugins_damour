BEGIN_PROVIDER [ double precision, OOOO, (spin_occ_num,spin_occ_num,spin_occ_num,spin_occ_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  OOOO(:,:,:,:) = dbERI(                                             &
      1:spin_occ_num,                                                &
      1:spin_occ_num,                                                &
      1:spin_occ_num,                                                &
      1:spin_occ_num   )
END_PROVIDER


BEGIN_PROVIDER [ double precision, OOOV, (spin_occ_num,spin_occ_num,spin_occ_num,spin_vir_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  OOOV(:,:,:,:) = dbERI(                                             &
      1:spin_occ_num,                                                &
      1:spin_occ_num,                                                &
      1:spin_occ_num,                                                &
      spin_occ_num+1:spin_mo_num)
END_PROVIDER


BEGIN_PROVIDER [ double precision, OOVO, (spin_occ_num,spin_occ_num,spin_vir_num,spin_occ_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  OOVO(:,:,:,:) = dbERI(                                             &
      1:spin_occ_num,                                                &
      1:spin_occ_num,                                                &
      spin_occ_num+1:spin_mo_num,                                    &
      1:spin_occ_num)                                                
END_PROVIDER


BEGIN_PROVIDER [ double precision, OVOO, (spin_occ_num,spin_vir_num,spin_occ_num,spin_occ_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  OVOO(:,:,:,:) = dbERI(                                             &
      1:spin_occ_num,                                                &
      spin_occ_num+1:spin_mo_num,                                    &
      1:spin_occ_num,                                                &
      1:spin_occ_num)                                                
END_PROVIDER


BEGIN_PROVIDER [ double precision, VOOO, (spin_vir_num,spin_occ_num,spin_occ_num,spin_occ_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  VOOO(:,:,:,:) = dbERI(                                             &
      spin_occ_num+1:spin_mo_num,                                    &
      1:spin_occ_num,                                                &
      1:spin_occ_num,                                                &
      1:spin_occ_num)                                                
END_PROVIDER


BEGIN_PROVIDER [ double precision, OOVV, (spin_occ_num,spin_occ_num,spin_vir_num,spin_vir_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  OOVV(:,:,:,:) = dbERI(                                             &
      1:spin_occ_num,                                                &
      1:spin_occ_num,                                                &
      spin_occ_num+1:spin_mo_num,                                    &
      spin_occ_num+1:spin_mo_num)
END_PROVIDER


BEGIN_PROVIDER [ double precision, OVOV, (spin_occ_num,spin_vir_num,spin_occ_num,spin_vir_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  OVOV(:,:,:,:) = dbERI(                                             &
      1:spin_occ_num,                                                &
      spin_occ_num+1:spin_mo_num,                                    &
      1:spin_occ_num,                                                &
      spin_occ_num+1:spin_mo_num)
END_PROVIDER


BEGIN_PROVIDER [ double precision, OVVO, (spin_occ_num,spin_vir_num,spin_vir_num,spin_occ_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  OVVO(:,:,:,:) = dbERI(                                             &
      1:spin_occ_num,                                                &
      spin_occ_num+1:spin_mo_num,                                    &
      spin_occ_num+1:spin_mo_num,                                    &
      1:spin_occ_num)
END_PROVIDER

BEGIN_PROVIDER [ double precision, OVVO_prime, (spin_occ_num,spin_vir_num,spin_vir_num,spin_occ_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  integer                        :: m,b,e,j
  do j=1,spin_occ_num
    do b=1,spin_vir_num
      do e=1,spin_vir_num
        do m=1,spin_occ_num
          OVVO_prime(m,e,b,j) = OVVO(m,b,e,j)
        enddo
      enddo
    enddo
  enddo
  FREE OVVO
END_PROVIDER

BEGIN_PROVIDER [ double precision, VOVO, (spin_vir_num,spin_occ_num,spin_vir_num,spin_occ_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  VOVO(:,:,:,:) = dbERI(                                             &
      spin_occ_num+1:spin_mo_num,                                    &
      1:spin_occ_num,                                                &
      spin_occ_num+1:spin_mo_num,                                    &
      1:spin_occ_num)
END_PROVIDER


BEGIN_PROVIDER [ double precision, VVOO, (spin_vir_num,spin_vir_num,spin_occ_num,spin_occ_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  VVOO(:,:,:,:) = dbERI(                                             &
      spin_occ_num+1:spin_mo_num,                                    &
      spin_occ_num+1:spin_mo_num,                                    &
      1:spin_occ_num,                                                &
      1:spin_occ_num)
END_PROVIDER


BEGIN_PROVIDER [ double precision, OVVV, (spin_occ_num,spin_vir_num,spin_vir_num,spin_vir_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  OVVV(:,:,:,:) = dbERI(                                             &
      1:spin_occ_num,                                                &
      spin_occ_num+1:spin_mo_num,                                    &
      spin_occ_num+1:spin_mo_num,                                    &
      spin_occ_num+1:spin_mo_num)
END_PROVIDER

BEGIN_PROVIDER [ double precision, VOVV, (spin_vir_num,spin_occ_num,spin_vir_num,spin_vir_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  VOVV(:,:,:,:) = dbERI(                                             &
      spin_occ_num+1:spin_mo_num,                                    &
      1:spin_occ_num,                                                &
      spin_occ_num+1:spin_mo_num,                                    &
      spin_occ_num+1:spin_mo_num)
END_PROVIDER


BEGIN_PROVIDER [ double precision, VVOV, (spin_vir_num,spin_vir_num,spin_occ_num,spin_vir_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  VVOV(:,:,:,:) = dbERI(                                             &
      spin_occ_num+1:spin_mo_num,                                    &
      spin_occ_num+1:spin_mo_num,                                    &
      1:spin_occ_num,                                                &
      spin_occ_num+1:spin_mo_num)
END_PROVIDER


BEGIN_PROVIDER [ double precision, VVVO, (spin_vir_num,spin_vir_num,spin_vir_num,spin_occ_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  VVVO(:,:,:,:) = dbERI(                                             &
      spin_occ_num+1:spin_mo_num,                                    &
      spin_occ_num+1:spin_mo_num,                                    &
      spin_occ_num+1:spin_mo_num,                                    &
      1:spin_occ_num)
END_PROVIDER

BEGIN_PROVIDER [ double precision, VVVV, (spin_vir_num,spin_vir_num,spin_vir_num,spin_vir_num) ]
  implicit none
  BEGIN_DOC
  END_DOC
  VVVV(:,:,:,:) = dbERI(                                             &
      spin_occ_num+1:spin_mo_num,                                    &
      spin_occ_num+1:spin_mo_num,                                    &
      spin_occ_num+1:spin_mo_num,                                    &
      spin_occ_num+1:spin_mo_num)
END_PROVIDER

