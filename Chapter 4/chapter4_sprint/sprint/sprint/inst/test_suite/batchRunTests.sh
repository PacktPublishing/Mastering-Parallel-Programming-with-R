#!/bin/bash
for (( c=1; c<=10; c++ ))
do
  /usr/local/bin/mpiexec -n 4 /usr/local/bin/R -f /Users/egrant1/Documents/SPRINT/workspace/sprint/trunk/test_suite/batch_ppam_once_test_suite.R  
done
