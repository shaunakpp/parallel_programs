#ifndef STRUCTS_H_
#define STRUCTS_H_

#define FILENAME "graphs/berlin52.tsp"
#define CITIES 52
#define TAG 100

// Default values
#define POPULATION 500
#define MUTATION_RATE 0.1
#define GENERATIONS  500
#define ELITE_CANDIDATES 100
#define MUTATION_SIZE 2
#define MIGRATIONS 500
#define CROSSOVERS 500
#define MIGRATION_SHARE 0.01

typedef struct{
  // population metrics
  int population_size;
  int generations;
  int elite_factor;
  // mutation factors
  double mutation_rate;
  int mutation_size;
  // migration factors
  int migrations;
  float migration_share;
  int crossovers;
}configuration;

typedef struct{
  double x;
  double y;
  unsigned int id;
}node;

typedef struct{
  int number_of_cities;
  node nodes[CITIES];
}city_graph;

typedef struct{
  int chromosomes[CITIES];
  int fitness;
  double probability;
  int index;
}candidate;

typedef struct{
  candidate *candidates;
  int number_of_candidates;
  double total_fitness;
}population;

configuration config = {POPULATION, GENERATIONS, ELITE_CANDIDATES, MUTATION_RATE, MUTATION_SIZE, MIGRATIONS, MIGRATION_SHARE, CROSSOVERS };
city_graph graph;

void init_population(population *pop);
void swap(int *a, int *b);
void calculate_fitness(population *pop);
int compare (const void *a, const void *b);
void sort(population *pop);
void crossover(population *pop);
void mutate(population *pop);
void migrate(population *pop, int rank, int numProcs);
void join_population(population *pop, int rank, int num_procs);
void print_candidate(candidate *cand);

#endif
