set terminal png size 1200,800 crop
set samples 10000
set output 'contour.jpg'

set xrange [0:99]
set yrange [0:99]
set size ratio 1
set tics out
set tics nomirror
set pm3d map

# palette
load 'blues.pal'


splot 'output.txt' using 1:2:5