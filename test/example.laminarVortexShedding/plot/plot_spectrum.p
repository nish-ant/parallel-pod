#!/usr/bin/gnuplot
# set term postscript eps size 1024, 720 color blacktext "Helvetica" 24
# set term epslatex color
set term eps

parentDIR="/scratch/nkumar001/rom4wt/POD/test"
datDIR=parentDIR."/example.laminarVortexShedding.run20220518"
podDIR=datDIR."/pod.m0.5.c8"
eigvalFile=podDIR."/chronos/chronos.bin" #- TODO: Read binary file

outDIR=parentDIR."/example.laminarVortexShedding.run20220518"
set output outDIR."/plot/plot_spectrum.eps"
set datafile commentschar '# '

set autoscale
unset log
unset label

set xtic auto
set ytic auto
unset title

set logscale y

set xlabel "Modes"
set ylabel "Eigenvalue"

set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 2 \
    pointtype 6 pointsize 0.5

plot eigvalFile every ::1::50 using ($0+1):1 notitle with linespoints linestyle 1