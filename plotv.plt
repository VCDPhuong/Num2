set terminal pdf
set output 'Data/initialvelo.pdf'
set xrange [-5:5]
binwidth = 2e-1
bin(x,width)=width*floor(x/width)
plot 'Data/Inivelocity.txt' using (bin($1,binwidth)):(1.0) smooth freq with lines,\

set output 'Data/Temperature.pdf'
unset xrange
set xlabel "Time (ps)"
plot 'Data/Argon.txt' u 1:2 w l t "Temperature",\
 94.4 w l

