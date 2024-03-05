#!/usr/bin/gnuplot
# set term postscript eps size 1024, 720 color blacktext "Helvetica" 24
set term eps

parentDIR="/scratch/nkumar001/rom4wt/POD/test/example.laminarVortexShedding.run20220518"
podDIR=parentDIR."/pod.m0.5.c8"
aFile=podDIR."/chronos/chronos.bin" #- TODO: Read binary file
tFile=parentDIR."/system/pod/snapshotTimes"

outDIR="/scratch/nkumar001/rom4wt/POD/test/example.laminarVortexShedding.run20220518"
set output outDIR."/plot/plot_apod.eps"
set datafile commentschar '# '

set autoscale
unset log
unset label

set xtic auto
set ytic auto
unset title

set xlabel "Time [s]"
set ylabel "a"

set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 2

# TODO: Plot based on input data
# plot for [i=1:2] '<paste '.tFile.' '.aFile.'chronos_'.i.'.dat' using 1:2 title 'a_'.i with lines linestyle i