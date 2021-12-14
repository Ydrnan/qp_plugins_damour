#!/bin/bash

source /home/yannd/Documents/Stage-M2/qp2/quantum_package.rc

echo "set contour" > test_functions/script_mc.gnu
echo "set cntrparam level incremental 0, 2, 1000" >> test_functions/script_mc.gnu
echo "set xrange [-1.5:4]" >> test_functions/script_mc.gnu
echo "set yrange [-3:4]" >> test_functions/script_mc.gnu
echo "set xlabel 'x'" >> test_functions/script_mc.gnu
echo "set ylabel 'y'" >> test_functions/script_mc.gnu
echo "set zlabel 'f(x,y)'" >> test_functions/script_mc.gnu
echo "splot sin(x+y)+(x-y)**2+(-1.5)*x+2.5*y+1 title 'McCormick function'" >> test_functions/script_mc.gnu

for i in {1..100}
do
		qp run test_mc > test_functions/data_mc${i}.dat > test_functions/opt_mc${i}.dat
		grep 'new pos' test_functions/opt_mc${i}.dat | awk '{printf "%-15s %-15s %-15s \n", $3, $4, $5}' > test_functions/data_mc${i}.dat
		echo "replot 'data_mc${i}.dat' u 1:2:3 title 'test $i' w linespoint pointsize 2" >> test_functions/script_mc.gnu
done
