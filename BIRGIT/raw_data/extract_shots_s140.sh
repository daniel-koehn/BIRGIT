#!/bin/bash

x1=4200
x2=4322

y1=1

step=1

x=$x1
y=$y1

while [ $x -lt $x2 ]
do

# extract shot
suwind < s140_geom_converted.su key=fldr min=$x max=$x > shots_s140/SH_s140_shot_$y.su

FILESIZE=$(stat -c%s "shots_s140/SH_s140_shot_$y.su")

echo "processing shot at $x "
#echo "file size =  $FILESIZE "

x=$[$x+$step]
y=$[$y+$step]

if [ $FILESIZE -lt $step ]; then
    echo "file size == 0 kB -> delete file and change counter"
    y=$[$y-$step]
    rm shots_s140/SH_s140_shot_$y.su
fi

done

