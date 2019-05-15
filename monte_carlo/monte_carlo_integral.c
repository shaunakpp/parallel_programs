#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

double function_to_integrate(double number){
  return(sqrt(number + sqrt(number)));
}

int main(int argc, char *argv[]){
  int i;
  int myid,numprocs,proc;
  MPI_Status status;
  int master = 0;
  long int iterations = 10000000;
  double local_sum = 0.0, sum;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  srand48(time(NULL)*myid);
  for(i = 0; i < iterations; i ++ ){
    local_sum += function_to_integrate((double)drand48());
  }
  local_sum = local_sum / iterations;
  if(myid == 0){
    sum = local_sum;
    for (proc=1; proc<numprocs; proc++) {
      MPI_Recv(&sum, 1, MPI_REAL, proc, 100, MPI_COMM_WORLD, &status);
      sum += local_sum;
    }
    sum = sum/numprocs;
    printf("Integration value: %f\n", sum);
  }else{
    MPI_Send(&local_sum, 1, MPI_REAL, master, 100, MPI_COMM_WORLD);
  }
  MPI_Finalize();
  return 0;
}
