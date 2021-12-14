#!/bin/bash

source /home/yannd/Documents/Stage-M2/qp2/quantum_package.rc

echo "set contour" > test_functions/script_gp.gnu
echo "set cntrparam level incremental 0, 100, 1000" >> test_functions/script_gp.gnu
echo "set xrange [-2:2]" >> test_functions/script_gp.gnu
echo "set yrange [-2:2]" >> test_functions/script_gp.gnu
echo "set xlabel 'x'" >> test_functions/script_gp.gnu
echo "set ylabel 'y'" >> test_functions/script_gp.gnu
echo "set zlabel 'f(x,y)'" >> test_functions/script_gp.gnu
echo "splot (1+(x+y+1)**2*(19-14*x+3*x**2+(-14)*y+6*x*y+3*y**2))*(30+(2*x-3*y)**2*(18-32*x+12*x**2+48*y+(-36)*x*y+27*y**2)) title 'Goldstein-Price function'" >> test_functions/script_gp.gnu

for i in {1..100}
do
		qp run test_gp > test_functions/data_gp${i}.dat > test_functions/opt_gp${i}.dat
		grep 'new pos' test_functions/opt_gp${i}.dat | awk '{printf "%-15s %-15s %-15s \n", $3, $4, $5}' > test_functions/data_gp${i}.dat
		echo "replot 'data_gp${i}.dat' u 1:2:3 title 'test $i' w linespoint pointsize 2" >> test_functions/script_gp.gnu
done
