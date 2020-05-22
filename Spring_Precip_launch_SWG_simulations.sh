#!/bin/sh -l
## File path needs to be adapted to local machine

start_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\n\nStarting script at: ${start_date}\n"


varname=rr
Lsim=122
daystart=01
nsim=1000
yy0=1950
alphaTN=0.5
alphacal=0.5
JOBID=precipAprJul
nn=5
fileanalo="data/analogues_slp_-50_30_30_70.txt"
##-----------------------------------------------------------------------

## Summer heatwaves
mostart=04
yy1=2018

staname=Orly

for alphacal in 0.5; do R CMD BATCH "--args ${varname} ${Lsim} ${mostart} ${daystart} ${nsim} ${yy0} ${yy1} ${fileanalo} ${JOBID} ${alphaTN} ${alphacal} ${nn}" SWG_ndaysCum.R log/${JOBID}_AnlmSa_${alphacal}_${alphaTN}_${nsim}_ndaysCum${nn}_${JOBID}.txt; done;

##-----------------------------------------------------------------------

start_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\n\nEnding script at: ${start_date}\n"
