#!/bin/bash

source /home/yannd/Documents/Stage-M2/qp2/quantum_package.rc

echo "set contour" > test_functions/script_poly.gnu
echo "set cntrparam level incremental 0, 100, 1000" >> test_functions/script_poly.gnu
echo "set xrange [-10:10]" >> test_functions/script_poly.gnu
echo "set yrange [-10:10]" >> test_functions/script_poly.gnu
echo "set xlabel 'x'" >> test_functions/script_poly.gnu
echo "set ylabel 'y'" >> test_functions/script_poly.gnu
echo "set zlabel 'f(x,y)'" >> test_functions/script_poly.gnu
echo "splot x**4 + y**4 title 'Polynomial function'" >> test_functions/script_poly.gnu

for i in {1..100}
do
		qp run test_poly > test_functions/data_poly${i}.dat > test_functions/opt_poly${i}.dat
		grep 'new pos' test_functions/opt_poly${i}.dat | awk '{printf "%-15s %-15s %-15s \n", $3, $4, $5}' > test_functions/data_poly${i}.dat
		echo "replot 'data_poly${i}.dat' u 1:2:3 title 'test $i' w linespoint pointsize 2" >> test_functions/script_poly.gnu
done
