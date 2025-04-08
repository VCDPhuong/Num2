clear
rm -r Data/
rm *.mod
mkdir Data/
gfortran Argon.f90 -fopenmp
./a.out
gnuplot plotE.plt
gnuplot plotv.plt