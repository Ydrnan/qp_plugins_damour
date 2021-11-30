#!/bin/bash

./auto_test.sh > result_test.dat
grep "ok" result_test.dat
echo ""
echo "======================="
echo "Number of tests: " $(grep -c "ok" result_test.dat) 
nb_error=$(grep -c "NOT ok" result_test.dat) 
nb_error2=$nb_error
echo "Number of failed tests: " $nb_error 
echo ""
