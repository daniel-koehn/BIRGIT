#!/bin/bash
DATA1="Kleinneudorf_true_FWI"

x1=1
x2=66          # no. of shots  (max(fldr)+1)
ntr=64         # no. of traces (max(tracf))
nt1=33333      # nt1 = number of samples _after_ interpolation (header value from DENISE synthetic data)
dt1=0.000030   # dt1 = sample interval from header of DENISE synthetic data (dt/10^6)

step=1
x=$x1

while [ $x -lt $x2 ]
do

# resample data 
# nt = number of samples _after_ interpolation (header value from DENISE synthetic data)
# dt = sample interval from header of DENISE synthetic data (dt/10^6)  
suresamp < tmp/DENISE_kleinneudorf_SH_y.su.shot$x nt=$nt1 dt=$dt1 > $DATA1/tmp1.su

# sort data
susort < $DATA1/tmp1.su tracf > $DATA1/DENISE_kleinneudorf_SH_y.su.shot$x


echo "processing shot ... $x"

x=$[$x+$step]
done

rm $DATA1/tmp*
