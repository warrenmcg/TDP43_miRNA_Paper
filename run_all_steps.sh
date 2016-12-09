#!/bin/bash

set -e

cd step_0/
bash step_0.sh &> log_step0.txt
cd ../step_1b/
bash step_1b.sh &> log_step1b.txt
cd ../step_2a/
bash step_2a.sh &> log_step2a.txt
cd ../step_2b/
bash step_2b.sh &> log_step2b.txt
cd ../step_2c/
bash step_2c.sh &> log_step2c.txt
cd ../step_3/
bash step_3.sh &> log_step3.txt
cd ../step_4/
bash step_4.sh &> log_step4.txt
cd ../step_5/
bash step_5.sh &> log_step5.txt
cd ../step_7/
bash step_7.sh &> log_step7.txt
cd ../step9/
bash step_9.sh &> log_step9.txt
cd ..

echo "Full run complete"
