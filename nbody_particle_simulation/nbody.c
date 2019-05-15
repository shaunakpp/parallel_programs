/* 2-D Gravitational N-Body simulation using Barnes-Hut algorithm
*  Shaunak Pagnis
*/
#include "mpi.h"
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "bhtree.c"

double min(double a, double b) {
    return a<b ? a : b;
}
double max(double a, double b) {
    return a>b ? a : b;
}
const double G = 6.673e-11;

void InitParticles(struct part *, int, int );


int main( int argc, char *argv[] )
{
  struct part   *particles;
  int           rank, size, npart, i;
  int           t_step;
  int           step;
  double        delta = 1.0;
  double        theta = 1.0;
  double        start_time;
  double        end_time;
  FILE *fp;
  char wfile[255];

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  srand48(time(NULL)*rank);


  npart = atoi(argv[1]) / size;
  t_step = atoi(argv[2]);

  particles = malloc(npart*sizeof(struct part));

  if(rank == 0){
    InitParticles(particles, npart, rank);
  }


  // if(rank == 0){
    sprintf(wfile,"timedat.%d", rank);
    fp=fopen(wfile, "w+");
  // }

  MPI_Bcast(&npart,1,MPI_INT, 0, MPI_COMM_WORLD);
  for(i=0; i < npart; i++){
    MPI_Bcast(&((particles)[i].mass), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&((particles)[i].x), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&((particles)[i].y), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&((particles)[i].vx), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&((particles)[i].vy), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&((particles)[i].fx), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&((particles)[i].fy), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  start_time = MPI_Wtime();
  for(step = 1; step<=t_step; step++){

    struct node * root;
    double xmin=0.0, xmax=0.0;
    double ymin=0.0, ymax=0.0;


    for(i=0; i<npart; i++){
      particles[i].fx=0.0;
      particles[i].fy=0.0;
      xmin=min(xmin,particles[i].x);
      xmax=max(xmax,particles[i].x);
      ymin=min(ymin,particles[i].y);
      ymax=max(ymax,particles[i].y);
    }

    root = create_node(particles+0,xmin,xmax,ymin,ymax);

    for(i=1; i<npart; i++){
      insert_particle(particles+i,root);
    }

    for(i=rank; i<npart; i+=size){
      calculate_forces(root,particles+i, G,theta);
    }

    for(i=0; i<npart; i++){
      MPI_Bcast(&(particles[i].fx),1,MPI_DOUBLE,i%size,MPI_COMM_WORLD);
      MPI_Bcast(&(particles[i].fy),1,MPI_DOUBLE,i%size,MPI_COMM_WORLD);
    }

    for(i=0; i<npart; i++){
      particles[i].x += delta*particles[i].vx;
      particles[i].y += delta*particles[i].vy;

      particles[i].vx += particles[i].fx * (delta/particles[i].mass);
      particles[i].vy += particles[i].fy * (delta/particles[i].mass);

    }
    destroy_tree(root);
    // if(rank == 0)
    {
      for (i=0; i<npart; i++)
      {
        fprintf(fp,"%d %d %d %f %f %f %f\n",rank, i, step, particles[i].x, particles[i].y,particles[i].vx, particles[i].vy);
      }
      fprintf(fp,"\n\n");
    // }
  }

  end_time=MPI_Wtime();
  if(rank==0){
    printf("Time required = %e seconds\n", end_time-start_time);
    fprintf(stderr, "%d,%d,%e\n", npart, size, end_time-start_time );
  }


  free(particles);
  MPI_Finalize();
  return 0;
} //end main

// Initialize particles at random positions with initial random velocity
void InitParticles(struct part *particles, int npart, int rank )
{
  int i;
  double mass = 5.0;
  for(i = 0; i < npart; i++)
  {
    particles[i].mass = mass;
    particles[i].vx = 0.5 - drand48();
    particles[i].vy = 0.5 - drand48();
    particles[i].x = 10 * (0.5 - drand48());
    particles[i].y = 10 * (0.5 - drand48());
    particles[i].fx = 0.0;
    particles[i].fy = 0.0;

    // printf("particle %d: x=%f y=%f, vx=%f vy=%f\n",i,particles[i].x,particles[i].y,particles[i].vx,particles[i].vy);
  }
}
