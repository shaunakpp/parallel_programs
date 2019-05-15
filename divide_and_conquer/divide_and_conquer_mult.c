#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

unsigned long long Power(unsigned long long my_mult, int my_rank, int num_procs, int exponent) {
  unsigned long long        mult;
  unsigned long long        temp;
  int        partner;
  unsigned   bitmask = 1;
  my_mult = pow(my_mult, exponent);
  mult = my_mult;
  while (bitmask < num_procs) {
    partner = my_rank ^ bitmask;
    // calculate local exponent value and add it to current local multiplication
    MPI_Sendrecv(&mult, 1, MPI_UNSIGNED_LONG_LONG, partner, 0,
      &temp, 1, MPI_UNSIGNED_LONG_LONG, partner, 0,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      mult *= temp;
      exponent *= 2;
      bitmask <<= 1;
    }
    return mult;
  }

int main(int argc, char* argv[]) {
  int num_procs, my_rank;
  unsigned long long my_mult;
  unsigned long long mult;
  int exponent = 2;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  exponent = 16/num_procs;
  if(my_rank == 0){
    my_mult = atoi(argv[1]);
    my_mult = (unsigned long long)my_mult;
  }
  // Send number to all procs
  MPI_Bcast(&my_mult, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  // Calculate the x ^ 16
  mult = Power(my_mult, my_rank, num_procs, exponent);
  if (my_rank == 0)
  printf("%llu ^ 16 = %llu\n", my_mult, mult);
  MPI_Finalize();
  return 0;
}
