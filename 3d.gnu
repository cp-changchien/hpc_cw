
set terminal png size 1200,800 crop
set output '3d.jpg'

set pm3d
set hidden3d

load 'blues.pal'


splot 'output.txt' using 1:2:5
