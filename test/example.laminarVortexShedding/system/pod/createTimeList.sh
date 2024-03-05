#!/bin/bash
#- File to generate custom time list

#- user input
tStart=150
tEnd=170
#- NOTE: nSkip= 0,1 --> No skipping
nSkip=0

#- Set locale 
#- See: https://stackoverflow.com/a/28238855/7473705
LC_NUMERIC=en_US.UTF-8

#- postProcessing path containing snapshots
# parentDIR=$(pwd)
# parentDIR="/scratch/nkumar001/OpenFOAM/nkumar001-2.4.0/run/LaminarVortexShedding.${ms}.run${dt}"
parentDIR=$1
snapsDir="${parentDIR}/postProcessing/internalField"
timeList="${parentDIR}/system/pod/snapshotTimes"

mkdir -p "${parentDIR}/system/pod/"

#- Generate list of times
/usr/bin/ls -A -1v $snapsDir | grep -E '^[0-9.]+$' | LC_ALL=C sort -g > $timeList 2>&1

#- Find lines corresponding to time range
#- See: https://stackoverflow.com/a/47541176/7473705
lineNumStart="$(grep -n -m 1 $tStart $timeList | cut -d: -f1)"
lineNumEnd="$(grep -n -m 1 $tEnd $timeList | cut -d: -f1)"

#- Trim to time range
#- See: https://stackoverflow.com/a/2237656/7473705
sed -i -n "${lineNumStart},${lineNumEnd}p" $timeList

#- Skip lines
#- See: https://superuser.com/a/396557/1140702
if [ "$nSkip" -gt 0 ]; then
    # sed -i -n "1p;0~${nSkip}p" $timeList
    sed -i -n "1~${nSkip}p" $timeList
fi