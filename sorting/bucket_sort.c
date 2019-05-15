// Homework 3: Parallel Bucket Sort
// Author: Shaunak Pagnis

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"


// Generates random numbers in the interval of total number of elements provided
void random_numbers_generator(long int array[], int no_of_elements, int numprocs){
  int i;
  int j;
  int group_size;
  group_size = no_of_elements/numprocs;

  for(j = 0; j < numprocs; j++){
    for(i = 0; i < group_size; i++) {
      array[j * group_size + i] = rand() % (2*group_size) + (j* 2 * group_size);
    }
  }
}

// comparison function used by stdlib quick sort function(qsort)
// this will be used for local sorting of the buckets
static int comparison(const void*x, const void*y){
  if ((*(long int *)x) > (*(long int *)y)){
    return (1);
  }
  if ((*(long int *)x) < (*(long int *)y)){
    return (-1);
  }
  return (0);
}


int main(int argc, char* argv[]){

  long int * array;
  long int * local;
  int my_id;
  int numprocs;
  int no_of_elements;
  int group_size;
  int i;
  int bucket_no;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);


  if(my_id == 0){
    no_of_elements = atoi(argv[1]);
    if(no_of_elements % numprocs != 0){
      printf("Please enter the number of elements which can be perfectly divided by number of processors\n");
      MPI_Abort(MPI_COMM_WORLD, 0);
      MPI_Finalize();
      return (1);
    }
    array = malloc(no_of_elements*sizeof(long int));
    random_numbers_generator(array, no_of_elements, numprocs);
    printf("The generated array is: \n");
    for(i=0;i<no_of_elements;i++){
      printf("%ld ", array[i]);
    }
    printf("\n");
  }

  // Send number of elements to all processors
  MPI_Bcast(&no_of_elements, 1, MPI_INT, 0, MPI_COMM_WORLD);
  group_size = no_of_elements / numprocs;

  local = malloc(group_size*sizeof(long int));
  // divide the array into n/p groups and send to all processors
  MPI_Scatter(array, group_size, MPI_INT, local, group_size, MPI_INT, 0, MPI_COMM_WORLD);

  // Send the key to the correct processor.
  if(my_id == 0){
    for(i=0;i < no_of_elements; i++){
      bucket_no = array[i]/(2* group_size);
      MPI_Send(&array[i], 1, MPI_LONG, bucket_no, 0, MPI_COMM_WORLD);
    }
  }
  // receive the correct key for local bucket
  for(i = 0; i < group_size; i++){
    MPI_Recv(&local[i],1, MPI_LONG, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status );
  }

  // sort the local buckets using qsort, which is a stdlib implementation of quick sort
  qsort((long int * ) local, group_size, sizeof(long int), comparison);

  // Gather the local groups into our global array and then print it
  MPI_Gather(local, group_size, MPI_LONG, array, group_size, MPI_LONG, 0, MPI_COMM_WORLD);

  if(my_id == 0){
    printf("The sorted array is: \n");
    for(i=0;i<no_of_elements;i++){
      printf("%ld ", array[i]);
    }
    printf("\n");
  }

  free(local);
  if (my_id==0) free(array);
  MPI_Finalize();
  return 0;
}
