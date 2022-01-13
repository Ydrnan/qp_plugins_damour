
FILE=
EZFIO1=${FILE}.ezfio

qp set_file $EZFIO1
qp set orbital optimization localization_method boys
qp set orbital optimization localization_max_nb_iter 1e4
qp set mo_class -d [1-] -a [-] -v []
qp run localization > $FILE.loc.out
qp set mo_class -c [1-] -a [-]

tar zcf $EZFIO1 ${EZFIO1}_save_loc_mos.tar.gz
