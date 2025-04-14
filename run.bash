clear
rm -r Data/
rm *.mod
rm ./a.out
mkdir Data/
gfortran Argon.f90 -fopenmp
./a.out
gnuplot plotE.plt