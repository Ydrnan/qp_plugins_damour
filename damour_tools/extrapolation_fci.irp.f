! One line command:
! echo "yourfile" | qp run extrapolation_fci 
! or
! qp run extrapolation_fci <<< yourfile

! Before doing the extrapolation, you must run:
! echo "yourfile.fci.out" | qp run extract_E_cipsi
! or
! qp run extract_E_cipsi <<< yorfilename.fci.out


program extract_E_cipsi

  implicit none

  character (len=100) :: qp_root
  character (len=100) :: filename

  write(*,*) 'file name ?'
  read(*,*) filename

  !path_hf = adjustl(adjustr(ezfio_dir)//'/hartree_fock/energy')
  !call execute_command_line('echo '//filename)
  
  ! on utilise une variable logique
  logical::existe
  
  ! on teste l'existence du fichier hartree_fock/energy
  inquire( file=filename, exist=existe)
 
  ! la réponse est dans la variable existe
  If    ( existe ) Then
     ! existe = .true.
     !write(*,*)'the file exists'
     !write(*,*)''
  ElseIf( .NOT. existe) Then
     ! existe = .false.
     write(*,*)'This file does not exist'
     STOP
  Endif

  CALL EXECUTE_COMMAND_LINE('python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/extrapolation_fci.py ' // filename)

end
