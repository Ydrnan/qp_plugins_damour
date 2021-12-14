#!/bin/bash

source /home/yannd/Documents/Stage-M2/qp2/quantum_package.rc

echo "set contour" > test_functions/script_booth.gnu
echo "set cntrparam level incremental 0, 100, 1000" >> test_functions/script_booth.gnu
echo "set xrange [-10:10]" >> test_functions/script_booth.gnu
echo "set yrange [-10:10]" >> test_functions/script_booth.gnu
echo "set xlabel 'x'" >> test_functions/script_booth.gnu
echo "set ylabel 'y'" >> test_functions/script_booth.gnu
echo "set zlabel 'f(x,y)'" >> test_functions/script_booth.gnu
echo "splot (x+2*y-7)**2+(2*x+y-5)**2 title 'Booth function'" >> test_functions/script_booth.gnu

for i in {1..100}
do
		qp run test_booth > test_functions/data_booth${i}.dat > test_functions/opt_booth${i}.dat
		grep 'new pos' test_functions/opt_booth${i}.dat | awk '{printf "%-15s %-15s %-15s \n", $3, $4, $5}' > test_functions/data_booth${i}.dat
		echo "replot 'data_booth${i}.dat' u 1:2:3 title 'test $i' w linespoint pointsize 2" >> test_functions/script_booth.gnu
done
