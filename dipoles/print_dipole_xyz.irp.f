program print_dipole_xyz
 implicit none
 read_wf = .True.
 TOUCH read_wf
 call print_dipole_moment_xyz_v2
 call print_transition_dipole_moment
 call print_oscillator_strength
 

end
