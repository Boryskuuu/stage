reset

set title "PROBA PRÉSENCE"
set xlabel "position( UNITÉÉÉÉÉ ???)"
set ylabel "proba"
set grid
set key autotitle columnhead

plot for [i=2:*] 'singlerun.txt' using 1:i with lines

