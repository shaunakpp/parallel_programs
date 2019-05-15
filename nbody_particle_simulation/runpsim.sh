#!/bin/bash
#Workflow script for particle sim
test_cases=2
while(($test_cases<=4))
do
  mydir=$PWD
  rm anim.gif
  rm -f timedat.*

  #execute code
  mpirun -np $test_cases ./nbodypipe $1 $2

  #performn visualization
  test=0;
  while (($test <= $2-1))
  do
  echo "set xrange [-20:20]
  set title \"$1 particles at timestep $test\"
  set yrange [-20:20]
  set grid
  set term gif size 640,480
  set output '$test.gif'
  plot \"timedat.0\" i $test u 4:5 pt 3  ps 1 t \"Node 0\";" >data_$test.gnu
  gnuplot data_$test.gnu
  let test=$test+1
  done


  #cleanup
  rm -f *.gnu
  ls *.gif | sort -nk1 | xargs ./gifmerge -10 -l0 >anim.vid
  rm -f *.gif
  mv anim.vid anim.gif
  let test_cases=$test_cases*2
done
