#!/bin/bash

source /home/yannd/Documents/Stage-M2/qp2/quantum_package.rc

echo "set contour" > test_functions/script_him.gnu
echo "set cntrparam level incremental 0, 10, 1000" >> test_functions/script_him.gnu
echo "set xrange [-5:5]" >> test_functions/script_him.gnu
echo "set yrange [-5:5]" >> test_functions/script_him.gnu
echo "set xlabel 'x'" >> test_functions/script_him.gnu
echo "set ylabel 'y'" >> test_functions/script_him.gnu
echo "set zlabel 'f(x,y)'" >> test_functions/script_him.gnu
echo "splot (x**2+y-11)**2+(x+y**2-7)**2 title 'Himmelblau function'" >> test_functions/script_him.gnu

for i in {1..100}
do
		qp run test_him > test_functions/data_him${i}.dat > test_functions/opt_him${i}.dat
		grep 'new pos' test_functions/opt_him${i}.dat | awk '{printf "%-15s %-15s %-15s \n", $3, $4, $5}' > test_functions/data_him${i}.dat
		echo "replot 'data_him${i}.dat' u 1:2:3 title 'test $i' w linespoint pointsize 2" >> test_functions/script_him.gnu
done
