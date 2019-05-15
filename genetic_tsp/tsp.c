#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <mpi.h>
#include "structs.h"
#include "graph.c"
#include "genesis.c"

int rank = 0;

int main(int argc, char **argv){
	population pop;
	int current_generation = 0;
	int num_procs = 0;
	double start, end;

	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	srand ( time(NULL)*rank );

	load_graph();
	init_population(&pop);

	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	while(config.generations > 0 && current_generation < config.generations){
		current_generation++;
		calculate_fitness(&pop);
		if(current_generation % 500 == 0){
			printf("Top candidate:\n");
			print_candidate(&pop.candidates[0]);
		}
		sort(&pop);
		crossover(&pop);
		mutate(&pop);
		if(current_generation % config.migrations == 0){
			migrate(&pop, rank, num_procs);
		}
	}

	join_population(&pop, rank, num_procs);
	sort(&pop);

	MPI_Barrier(MPI_COMM_WORLD);
	end = MPI_Wtime();

	if(rank == 0) {
		printf("Execution time required: %f\n", end - start);
	}

	free(pop.candidates);
	MPI_Finalize();
	return 0;
}
