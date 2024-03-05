#!/usr/bin/gnuplot
# set term postscript eps size 1024, 720 color blacktext "Helvetica" 24
set term eps


parentDIR="/scratch/nkumar001/rom4wt/POD/test/example.laminarVortexShedding.run20220518"
datFile=parentDIR."/postProcessing/forceCoeffs/0/forceCoeffs.dat"

outDIR="/scratch/nkumar001/rom4wt/POD/test/example.laminarVortexShedding.run20220518"
set output outDIR."/plot/plot_ClCd_Re100.eps"
set datafile commentschar '# '

set autoscale
unset log
unset label

set xtic auto
set ytic auto
unset title

set xlabel "Time [s]"
set ylabel "C_D, C_L"

set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 2
set style line 2 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 2

plot datFile using 1:3 title "C_D" with lines linestyle 1, \
     ''      using 1:4 title "C_L" with lines linestyle 2