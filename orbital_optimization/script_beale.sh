#!/bin/bash

source /home/yannd/Documents/Stage-M2/qp2/quantum_package.rc

echo "set contour" > test_functions/script_beale.gnu
echo "set cntrparam level incremental 0, 10, 1000" >> test_functions/script_beale.gnu
echo "set xrange [-4.5:4.5]" >> test_functions/script_beale.gnu
echo "set yrange [-4.5:4.5]" >> test_functions/script_beale.gnu
echo "set xlabel 'x'" >> test_functions/script_beale.gnu
echo "set ylabel 'y'" >> test_functions/script_beale.gnu
echo "set zlabel 'f(x,y)'" >> test_functions/script_beale.gnu
echo "splot (1.5-x+x*y)**2+(2.25-x+x*y**2)**2+(2.625-x+x*y**3)**2 title 'Beale function'" >> test_functions/script_beale.gnu

for i in {1..100}
do
		qp run test_beale > test_functions/data_beale${i}.dat > test_functions/opt_beale${i}.dat
		grep 'new pos' test_functions/opt_beale${i}.dat | awk '{printf "%-15s %-15s %-15s \n", $3, $4, $5}' > test_functions/data_beale${i}.dat
		echo "replot 'data_beale${i}.dat' u 1:2:3 title 'test $i' w linespoint pointsize 2" >> test_functions/script_beale.gnu
done
