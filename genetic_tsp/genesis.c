#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include "structs.h"

// Initialize population with random chromosome values
void init_population(population *pop){

	pop->number_of_candidates = config.population_size;
	pop->candidates = (candidate*)malloc(pop->number_of_candidates * sizeof(candidate));

	int i, j, t;
	unsigned int k;
	for(i = 0; i < pop->number_of_candidates; i++){
		for(j = 0; j < CITIES; j++){
			pop->candidates[i].chromosomes[j] = j;
		}
		// randomize the tour by shuffling values
		srand48(time(NULL)*i);
		for (j = CITIES - 1; j > 0; j--) {
			k = (unsigned int) (drand48()*(j+1));
			t = pop->candidates[i].chromosomes[k];
			pop->candidates[i].chromosomes[k] = pop->candidates[i].chromosomes[j];
			pop->candidates[i].chromosomes[j] = t;
		}
	}
	return;
}

// Calculate fitness of each candidate based on Euclidean distance
void calculate_fitness(population *pop){
	int i,j;
	int total_fitness = 0.0;
	int fitness = 0.0;
	int city, next_city;
	int xd, yd;
	for(i=0; i<pop->number_of_candidates; i++){
		fitness = 0.0;

		for (j = 0; j < CITIES-1; j++){
			city = pop->candidates[i].chromosomes[j];
			next_city = pop->candidates[i].chromosomes[j+1];
			xd = graph.nodes[city].x - graph.nodes[next_city].x;
			yd = graph.nodes[city].y - graph.nodes[next_city].y;
			fitness = fitness +  rint(sqrt(xd*xd + yd*yd));
		}
		pop->candidates[i].fitness = fitness;
		total_fitness += fitness;
	}
	pop->total_fitness = total_fitness;
	return;
}

void swap (int *a, int *b){
	int temp = *a;
	*a = *b;
	*b = temp;
}

int compare (const void *a, const void *b)
{
	candidate * population_1 = (candidate *)a;
	candidate * population_2 = (candidate *)b;
	return ( population_1->fitness - population_2->fitness );
}

//mutate the population based on the elitism factor
void mutate(population *pop){
	candidate *cand = NULL;
	int swapped_list[CITIES];
	int mutation_rate = config.mutation_rate * 100;
	int index_1, index_2;
	int i, j;

	for(i = config.elite_factor; i < pop->number_of_candidates; i++){
		if(rand()%100 > mutation_rate){
			continue;
		}

		memset(swapped_list, 0, sizeof(swapped_list));
		cand = &pop->candidates[i];

		for(j=0; j<config.mutation_size; j++){
			index_1 = rand() % CITIES;

			while(swapped_list[index_1] !=0){
				if(index_1 + 1 < CITIES){
					index_1++;
				}
				else{
					index_1 = 0;
				}
			}
			swapped_list[index_1] = 1;

			index_2 = rand() % CITIES;

			while(swapped_list[index_2] !=0){
				if(index_2 + 1 < CITIES){
					index_2++;
				}
				else{
					index_2 = 0;
				}
			}
			swapped_list[index_2] = 1;
			swap(&cand->chromosomes[index_1], &cand->chromosomes[index_2]);
		}
	}
	return;
}

//sort the population by increasing fitness
void sort(population *pop){
	qsort(pop->candidates, pop->number_of_candidates, sizeof(candidate), compare);
	return;
}

// Using the PMX crossover, generate new children
void crossover(population *pop){

	int i, j;
	double random_1, random_2 ;
	int cross_pone, cross_ptwo;
	double offset_one = 0.0;
	double offset_two = 0.0;
	int pick_one = 0, pick_two = 0;
	int count;
	int max = CITIES -1;
	int min = 1;


	for (i=0; i < config.population_size; i++){
		pop->candidates[i].probability = pop->candidates[i].fitness/pop->total_fitness;
	}

	for (count = 0 ; count < config.crossovers; count= count + 2){
		random_1 = rand() / (double) RAND_MAX;
		random_2 = rand() / (double) RAND_MAX;
		for (i = 0; i < config.population_size; i++) {
			offset_one += pop->candidates[i].probability;
			if (random_1 < offset_one) {
				pick_one = i;
				break;
			}
		}
		for (i = 0; i < config.population_size; i++) {
			offset_two += pop->candidates[i].probability;
			if (random_2 < offset_two) {
				pick_two = i;
				break;
			}
		}

		cross_pone = (max - min + 1)*(double)rand()/RAND_MAX + min;
		cross_ptwo = (max - min + 1)*(double)rand()/RAND_MAX + min;

		if(cross_ptwo < cross_pone){
			swap(&cross_pone, &cross_ptwo);
		}
		int child_1[CITIES] = {0};
		int child_2[CITIES] = {0};

		for (i=0; i<CITIES; i++){
			child_1[i] = pop->candidates[pick_one].chromosomes[i];
			child_2[i] = pop->candidates[pick_two].chromosomes[i];
		}
		int *vec_1 = (int*)malloc((cross_ptwo-cross_pone)*sizeof(int));
		int *vec_2 = (int*)malloc((cross_ptwo-cross_pone)*sizeof(int));
		int num=0;

		for (i = cross_pone; i < cross_ptwo; i++ ){
			swap(&child_1[i], &child_2[i]);
			vec_1[num] = child_1[i];
			vec_2[num] = child_2[i];
			num++;
		}

		for(i = 0; i < cross_pone; i++){
			child_1[i] = pop->candidates[pick_one].chromosomes[i];
			child_2[i] = pop->candidates[pick_two].chromosomes[i];
		}
		for(i = cross_ptwo; i<CITIES; i++){
			child_1[i] = pop->candidates[pick_one].chromosomes[i];
			child_2[i] = pop->candidates[pick_two].chromosomes[i];
		}

		int flag=0;
		for (i = 0; i < num; i++ ){
			if(flag==1){
				i=i-1;
				flag=0;
			}
			for (j = 0; j < num; j++ ){
				if((vec_1[i]==vec_2[j])&&(vec_1[i]!=vec_2[i] )){
					swap(&vec_1[i], &vec_1[j]);
					flag=1;
					continue;
				}
			}
		}
		for (i=0; i<num; i++){
			for(j = 0; j<cross_pone; j++){
				if(vec_1[i]==child_1[j]){
					child_1[j] = vec_2[i];
				}
				if(vec_2[i]==child_2[j]){
					child_2[j] = vec_1[i];
				}

			}
			for(j = cross_ptwo; j<CITIES; j++){
				if(vec_1[i]==child_1[j]){
					child_1[j] = vec_2[i];
				}
				if(vec_2[i]==child_2[j]){
					child_2[j] = vec_1[i];
				}
			}
		}
		for (i=0; i<CITIES; i++){
			pop->candidates[config.population_size-count-1].chromosomes[i] = child_1[i];
			pop->candidates[config.population_size-count-2].chromosomes[i] = child_2[i];
		}
		free(vec_1);
		free(vec_2);
	}
	return;
}

// Migrate population across processors in a Ring baed topology
void migrate(population *pop, int rank, int numProcs){

	candidate *outgoing1 = NULL, *outgoing2 = NULL;
	candidate *incoming1 = NULL, *incoming2 = NULL;
	int i=0, j=-1;
	MPI_Status status;
	int number_of_migrants = floor(pop->number_of_candidates * config.migration_share);
	int size_of_migrants = sizeof(candidate) * number_of_migrants;

	i = rand() % (pop->number_of_candidates - number_of_migrants);
	outgoing1 = &pop->candidates[i];
	incoming1 = (candidate *)malloc(size_of_migrants);

	if(rank > 0 && rank < numProcs-1){
		do {
			j = rand() % (pop->number_of_candidates - number_of_migrants);
		} while(j == i);

		outgoing2 = &pop->candidates[j];
		incoming2 = (candidate *)malloc(size_of_migrants);
	}

	//Build a ring of processors and exchange migrant in clock wise direction
	if (rank % 2 != 0){
		if(rank < numProcs - 1){
			MPI_Send(outgoing1, size_of_migrants, MPI_CHAR, rank + 1, TAG, MPI_COMM_WORLD);
			MPI_Recv(incoming1, size_of_migrants, MPI_CHAR, rank + 1, TAG, MPI_COMM_WORLD, &status );
		}

		if(rank > 0){
			MPI_Send((rank < numProcs - 1? outgoing2 : outgoing1), size_of_migrants, MPI_CHAR, rank-1, TAG, MPI_COMM_WORLD);
			MPI_Recv((rank < numProcs - 1? incoming2 : incoming1), size_of_migrants, MPI_CHAR, rank-1, TAG, MPI_COMM_WORLD, &status );
		}

	} else {
		if(rank > 0){
			MPI_Recv(incoming1, size_of_migrants, MPI_CHAR, rank-1, TAG, MPI_COMM_WORLD, &status );
			MPI_Send(outgoing1, size_of_migrants, MPI_CHAR, rank-1, TAG, MPI_COMM_WORLD);
		}

		if(rank < numProcs-1){
			MPI_Recv((rank > 0? incoming2 : incoming1), size_of_migrants, MPI_CHAR, rank+1, TAG, MPI_COMM_WORLD, &status );
			MPI_Send((rank > 0? outgoing2 : outgoing1), size_of_migrants, MPI_CHAR, rank+1, TAG, MPI_COMM_WORLD);
		}
	}

	// move imigrants to our current population
	if(numProcs > 1){
		memcpy(&pop->candidates[i], incoming1, size_of_migrants);
	}

	if(rank > 0 && rank < numProcs-1){
		memcpy(&pop->candidates[j], incoming2, size_of_migrants);
	}

	free(incoming1);
	free(incoming2);
	return;
}

// Join the best candidates to create unified population
void join_population(population *pop, int rank, int num_procs){
	MPI_Status status;
	int i, j;
	int candidates_per_proc = floor(pop->number_of_candidates / num_procs);
	int k = candidates_per_proc;

	if(rank == 0){
		// get the top candidates from other populations
		for(i=1; i<num_procs; i++){
			for(j=0; j<candidates_per_proc; j++){
				MPI_Recv(&pop->candidates[k], sizeof(candidate), MPI_CHAR, i, TAG, MPI_COMM_WORLD, &status);
				k++;
			}
		}
	}else{
		for(j=0; j<candidates_per_proc; j++){
			MPI_Send(&pop->candidates[j], sizeof(candidate), MPI_CHAR, 0, TAG, MPI_COMM_WORLD);
			k++;
		}
	}
	return;
}

// Used for debugging
void print_candidate(candidate *cand){
	int i;
	printf("candidate_fitness: %d\n", cand->fitness);
	printf("chromosomes:\n");
	for(i=0; i<CITIES-1; i++){
		printf("<%d>->", cand->chromosomes[i]);
	}
	printf("<%d>\n",cand->chromosomes[CITIES-1]);
	printf("\n\n");
	return;
}
