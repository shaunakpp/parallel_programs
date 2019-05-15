// Homework 3: Parallel Sample Sort
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
  int i, j ,k, c;
  int bucket_no;
  MPI_Status status;
  int * local_splitter;
  int * global_splitter;
  int * global_buckets;
  int * local_bucket;
  int * bucket_size;
  int * op_size;
  int * op;

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
      /* printf("%ld ", array[i]); */
    }
    /* printf("\n"); */
  }

  // Send number of elements to all processors
  MPI_Bcast(&no_of_elements, 1, MPI_INT, 0, MPI_COMM_WORLD);
  group_size = no_of_elements / numprocs;

  local = malloc(group_size*sizeof(long int));
  // divide the array into n/p groups and send to all processors
  MPI_Scatter(array, group_size, MPI_INT, local, group_size, MPI_INT, 0, MPI_COMM_WORLD);

  // sort the local buckets using qsort, which is a stdlib implementation of quick sort
  qsort((long int * ) local, group_size, sizeof(long int), comparison);

  local_splitter = (int *) malloc(sizeof(long int) * numprocs-1);
  for(i=0; i < numprocs-1; i++){
    local_splitter[i] = local[no_of_elements/(numprocs * numprocs) * (i + 1)];
  }

  global_splitter = (int *) malloc(sizeof(long int) * numprocs * (numprocs - 1));
  MPI_Gather(local_splitter, numprocs-1, MPI_LONG, global_splitter, numprocs-1, MPI_LONG, 0, MPI_COMM_WORLD);

  if(my_id == 0){
    qsort((char*) global_splitter, numprocs * (numprocs - 1), sizeof(long int), comparison);
    for(i=0;i<numprocs-1;i++){
      local_splitter[i] = global_splitter[(numprocs - 1) * (i + 1)];
    }
  }

  MPI_Bcast(local_splitter, numprocs-1, MPI_LONG, 0, MPI_COMM_WORLD);

  global_buckets = (int *) malloc (sizeof (long int) * (no_of_elements + numprocs));

  j = 0;
  k = 1;
  for(i = 0; i < group_size; i++){
    if(j < (numprocs - 1)){
      if(local[i] < local_splitter[j]){
        global_buckets[((group_size + 1) * j) + k++] = local[i];
      }else{
        global_buckets[((group_size + 1) * j)] = k-1;
        k = 1;
        j++;
        i--;
      }
    }else{
      global_buckets[((group_size + 1) * j) + k++] = local[i];
    }
  }
  global_buckets[((group_size + 1) * j)] = k - 1;
  bucket_size = (int *) malloc(sizeof (long int) * (no_of_elements + numprocs));

  MPI_Alltoall(global_buckets, group_size + 1, MPI_LONG, bucket_size, group_size + 1, MPI_LONG, MPI_COMM_WORLD);

  local_bucket = (int *) malloc(sizeof(long int) * 2 * group_size);

  c = 1;
  for(j = 0; j < numprocs; j++){
    k = 1;
    for(i = 0; i < bucket_size[(group_size + 1) * j]; i++){
      local_bucket[c++] = bucket_size[((group_size + 1) * j) + k++];
    }
  }
  local_bucket[0] = c-1;

  qsort((char *) &local_bucket[1], c-1, sizeof(long int), comparison);
  if(my_id == 0) {
    op_size = malloc (sizeof(long int) * 2 * no_of_elements);
    op = malloc (sizeof (long int) * no_of_elements);
  }

  // Gather the local groups into our global array and then print it
  MPI_Gather(local_bucket, 2*group_size, MPI_LONG, op_size, 2*group_size, MPI_LONG, 0, MPI_COMM_WORLD);

  if(my_id == 0){
    c = 0;
    for(j=0; j<numprocs; j++){
      k = 1;
      for(i=0; i<op_size[(group_size) * j]; i++) {
        op[c++] = op_size[(group_size) * j + k++];
      }
    }
    free(array);
    free(op_size);
    free(op);
  }

  free(local);
  free(local_splitter);
  free(global_splitter);
  free(local_bucket);
  free(global_buckets);
  free(bucket_size);

  MPI_Finalize();
  return 0;
}
