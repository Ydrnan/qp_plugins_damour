program print_dipole_xyz

 implicit none

 BEGIN_DOC
 ! To print:
 ! - Dipole moments of the different states
 ! - Transition dipole moments
 ! - Oscillator strengths
 END_DOC

 read_wf = .True.
 TOUCH read_wf

 call print_dipole_moment
 call print_transition_dipole_moment
 call print_oscillator_strength

end
