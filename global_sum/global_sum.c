#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include <string.h>

#define MAXSIZE 1000

int main(int argc, char **argv)
{
	int myid, numprocs;
	int data[MAXSIZE];
	int i, x, high, offset, partial_sum, partial_max, mymax = 0, max = 0,  myresult = 0, result = 0;
	char fn[255];
	FILE *fp;
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	x = MAXSIZE/numprocs; /* Add my portion of data */

	if (myid == 0) {

		/* Open input file and initialize data */
		strcpy(fn, "rand_data.txt");
		if ((fp = fopen(fn,"r")) == NULL) {
			printf("Cannot open the input file: %s\n\n", fn);
			exit(1);
		}
		offset = x;
		for(i = 0; i < MAXSIZE; i++) fscanf(fp,"%d", &data[i]);

		// Divide data into chunk and send it to processors
		for(i = 1; i < numprocs; i++){
			MPI_Send(&offset,1 , MPI_INT, i, 1, MPI_COMM_WORLD);
			if(i == numprocs - 1){
				MPI_Send(&data[offset], (x + (MAXSIZE % numprocs)), MPI_INT, i, 0, MPI_COMM_WORLD);
			}else {
				MPI_Send(&data[offset], x, MPI_INT, i, 0, MPI_COMM_WORLD);
				offset += x;
			}
		}

		// Calculate sum of first chunk
		for(i = 0; i < x; i ++){
			result += data[i];
		}

		// Calculate Max value from first chunk
		for(i = 0; i < x; i ++){
			if(data[i] > max){
				max = data[i];
			}
		}

		// Receive partial_sums from other proccesors
		for(i = 1; i < numprocs; i++){
			MPI_Recv(&partial_sum, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
			result += partial_sum;
		}

		// Receive partial_sums from other proccesors
		for(i = 1; i < numprocs; i++){
			MPI_Recv(&partial_max, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
			if(partial_max > max){
				max = partial_max;
			}
		}
	} else {

		// Receive offset value
		MPI_Recv(&offset, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

		// Calculate limit for iterating over data
		if(offset == (MAXSIZE - x - (MAXSIZE % numprocs))){
			high = MAXSIZE;
		}else{
			high = offset + x;
		}

		// Receive data
		if(myid == numprocs - 1){
			MPI_Recv(&data[offset], (x + (MAXSIZE % numprocs)), MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		}else {
			MPI_Recv(&data[offset], x, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		}

		// Calculate partial sum
		for(i = offset; i < high; i++){
			myresult += data[i];
		}

		for(i = offset; i < high; i ++){
			if(data[i] > mymax){
				mymax = data[i];
			}
		}
		// send partial max result back to proccesor 1
		MPI_Send(&mymax, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		// send partial sum result back to proccesor 0
		MPI_Send(&myresult, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	if (myid == 0) printf("The sum is %d.\n", result);
	if (myid == 0) printf("The max number is %d.\n", max);
	MPI_Finalize();
	return 0;
}
