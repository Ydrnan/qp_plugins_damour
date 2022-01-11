source /home/yannd/Documents/Stage-M2/qp2/quantum_package.rc

qp set coupled_cluster sccd_method bi_int
qp run spin_orb_ccd > toto_bi_int

qp set coupled_cluster sccd_method guess_mp2
qp run spin_orb_ccd > toto_guess_mp2

qp set coupled_cluster sccd_method estimated_e
qp run spin_orb_ccd > toto_estimated_e

grep "Result" toto_*

