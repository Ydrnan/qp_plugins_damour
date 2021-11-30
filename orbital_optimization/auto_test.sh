#!/bin/bash

#RED='\033[1;31m'
#GREEN='\033[1;32m'
#NC='\033[0m' # No Color

#	source $QP_ROOT/quantum_package.rc
source ~/Documents/Stage-M2/qp2/quantum_package.rc

# script to check the results and compare them with the result of reference

#touch h.ezfio
#rm -r *.ezfio
#
rm -r reference_results
tar zxf reference_results.tar.gz
cd reference_results
cp ../methods_check.dat .
cp ../diff.py .

list_save_ezfio=$(ls -d *save.ezfio)

for file in ${list_save_ezfio}
do
    echo $file
	test0=$(sed 's/save/test/' <<< ${file})	
	echo $test0
	cp -r ${file} ${test0}
done

#list_ezfio=$(ls -d *test.ezfio)
ezfios=$(cat methods_check.dat | awk '{print $1}')
list_ezfio=($ezfios)
methods=$(cat methods_check.dat | awk '{print $2}')
list_method=($methods)
keyword=$(cat methods_check.dat | awk '{print $3}')
list_keyword=($keyword)
res=$(cat methods_check.dat | awk '{print $4}')
list_res=($res)

for ((i=0; i<${#list_method[@]}; i++))
do
	list_keyword_space[$i]=$(sed 's/_/ /g' <<< ${list_keyword[$i]})
	#echo ${list_keyword[$i]}

	#echo "toto" ${list_keyword_space[$i]}
done

not_ok=0

for ((i=0; i<${#list_ezfio[@]}; i++))
do
    echo "###############################"
    echo ${list_ezfio[$i]}
	qp set_file ${list_ezfio[$i]}

    save=$(sed "s/test/save/" <<< ${list_ezfio[$i]})
    rm -r ${list_ezfio[$i]}
    cp -r $save ${list_ezfio[$i]}

    #echo $i
    echo ${list_method[$i]}
    #echo "key"  ${list_keyword_space[$i]}
    v1=$(qp run ${list_method[$i]} | grep "${list_keyword_space[$i]}" | tail -n 1 | awk '{print $NF}')
    #grep "${list_method[$i]}" methods_check.dat  | awk '{print $NF}' > tmp_v2  #awk -vi=${list_columns} '{print $i}')
    v2=${list_res[$i]}
    #v1=$(cat tmp_v1)
    #v2=$(cat tmp_v2)
    echo "v1" $v1
    echo "v2" $v2
    ## print
    python3 diff.py $v1 $v2
    is_ok=$(python3 diff.py $v1 $v2 | grep "Valid" | awk '{print $2}')
    if [ "$is_ok" = "True" ]
    then
            echo "${list_ezfio[$i]} ${list_method[$i]} ok"
    else
            echo "${list_ezfio[$i]} ${list_method[$i]} NOT ok"
			not_ok=$((not_ok+1))
	fi

    echo ""
	echo "Number of errors: " $not_ok

done

#for ezfio in $list_ezfio
#do
#    echo "###############################"
#    echo $ezfio
#	qp set_file $ezfio
# 
#	for ((i=0; i<${#list_method[@]}; i++))
#	do
#		save=$(sed "s/test/save/" <<< $ezfio)
#		rm -r $ezfio
#		cp -r $save $ezfio
#
#		#echo $i
#		echo ${list_method[$i]}
#		#echo "key"  ${list_keyword_space[$i]}
#		v1=$(qp run ${list_method[$i]}| grep "${list_keyword_space[$i]}" | tail -n 1 | awk '{print $NF}')
#		#grep "${list_method[$i]}" methods_check.dat  | awk '{print $NF}' > tmp_v2  #awk -vi=${list_columns} '{print $i}')
#		v2=${list_res[$i]}
#		#v1=$(cat tmp_v1)
#		#v2=$(cat tmp_v2)
#		echo "v1" $v1
#		echo "v2" $v2
#		## print
#        python3 diff.py $v1 $v2
#        is_ok=$(python3 diff.py $v1 $v2 | grep "Valid" | awk '{print $2}')
#        if [ "$is_ok" = "True" ]
#        then
#                echo "$ezfio ${list_method[$i]} ok"
#        else
#                echo "$ezfio ${list_method[$i]} NOT ok"
#        fi
#
#		echo ""
#		echo "===================================="
#		echo ""
#
#	done 
#	echo ""
#done

#for ezfio in $list_ezfio
#do 
#
#    qp set_file $ezfio
#	echo $ezfio
#	echo "scf check"
#	
#	v1=$(qp run scf | grep "SCF energy"  )
#	v2=$(grep "$ezfio SCF energy" | awk '{print $4}')
#	# print
#	python3 diff.py $v1 $v2
#	is_ok=$(python3 diff.py $v1 $v2 | grep "Valid" | awk '{print $2}')
#	if [ "$is_ok" = "True" ]
#	then
#			echo "$ezfio SCF ok"
#	else
#			echo "$ezfio SCF NOT ok"
#	fi
#done
