set terminal pdf
set output 'Data/Energy.pdf'
set xlabel "Time"
plot  "Data/Energy.txt" u 1:2 w l title "Kinetic Energy",\
      "Data/Energy.txt" u 1:3 w l title "Potential energy",\
      "Data/Energy.txt" u 1:4 w l title "Total Energy",\
      1.380649*10.**(-16.)*94.4*3./2.*864*1e12 title "T = 94.4 K"

set output "Data/Temp.pdf"
plot "Data/Energy.txt" u 1:5 w l,\
      94.4 title "T = 94.4 K"
set xlabel "Time [ps]"
set ylabel "Temperature