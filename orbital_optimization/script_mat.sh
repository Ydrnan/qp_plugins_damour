#!/bin/bash

source /home/yannd/Documents/Stage-M2/qp2/quantum_package.rc

echo "set contour" > test_functions/script_mat.gnu
echo "set cntrparam level incremental 0, 10, 100" >> test_functions/script_mat.gnu
echo "set xrange [-10:10]" >> test_functions/script_mat.gnu
echo "set yrange [-10:10]" >> test_functions/script_mat.gnu
echo "set xlabel 'x'" >> test_functions/script_mat.gnu
echo "set ylabel 'y'" >> test_functions/script_mat.gnu
echo "set zlabel 'f(x,y)'" >> test_functions/script_mat.gnu
echo "splot 0.26*(x**2+y**2)-0.48*x*y title 'Matyas function'" >> test_functions/script_mat.gnu

for i in {1..100}
do
		qp run test_mat > test_functions/data_mat${i}.dat > test_functions/opt_mat${i}.dat
		grep 'new pos' test_functions/opt_mat${i}.dat | awk '{printf "%-15s %-15s %-15s \n", $3, $4, $5}' > test_functions/data_mat${i}.dat
		echo "replot 'data_mat${i}.dat' u 1:2:3 title 'test $i' w linespoint pointsize 2" >> test_functions/script_mat.gnu
done
