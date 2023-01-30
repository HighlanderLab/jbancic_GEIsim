#!/bin/bash

# no spaces when giving commands
SCENARIOS='02RUNMEwithGEI_scenarios' 

# run scenarios 
for i in $SCENARIOS;
do
  	qsub -t 1-10 ${i}.R
done

echo Jobs submitted
