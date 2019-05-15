#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "particle.h"
#include "bhtree.h"

enum quadrant{NE,NW,SW,SE};
enum quadrant get_quadrant(double x, double y, double xmin, double xmax, double ymin, double ymax);

//functions for tree manipulations
struct node *create_node(struct part * particle, double xmin, double xmax, double ymin, double ymax);
void insert_particle(struct part * particle, struct node * bh_node);
void destroy_tree(struct node * bh_node);

// functions for force calculations
void update_mass_of_center(struct node * bh_node, struct part * particle);
double update_center(double center, double position, double total_mass, double mass);
void calculate_forces(struct node * bh_node, struct part * particle, double G, double thetha );

// Returns the box for a particle based on the position
enum quadrant get_quadrant(double x, double y, double xmin, double xmax, double ymin, double ymax){
    double midx, midy;

    midx = xmin + 0.5*(xmax-xmin);
    midy = ymin + 0.5*(ymax-ymin);

    if(y>midy){
        if(x>midx)
            return NE;
        else
            return NW;
    }
    else{
        if(x>midx)
            return SE;
        else
            return SW;
    }

}

// Create root node for the tree and add the first particle
struct node *create_node(struct part *particle, double xmin, double xmax, double ymin, double ymax){
    struct node *root;
    if(!(root=malloc(sizeof(struct node)))){
        printf("Can't allocate node\n");
    }

    root->total_mass = particle->mass;
    root->centerx = particle->x;
    root->centery = particle->y;
    root->xmin = xmin;
    root->xmax = xmax;
    root->ymin = ymin;
    root->ymax = ymax;
    root->box_size = sqrt(( pow(xmax - xmin, 2) + pow(ymax - ymin, 2) ));
    root->particle = particle;
    root->NE = NULL;
    root->NW = NULL;
    root->SW = NULL;
    root->SE = NULL;

    return root;
}

// Insert the particle into correct position in the tree
void insert_particle(struct part *new_particle, struct node *bh_node){

    enum quadrant existing_quadrant, new_quadrant;
    double xmid, ymid;
    xmid = bh_node->xmin + 0.5*(bh_node->xmax - bh_node->xmin);
    ymid = bh_node->ymin + 0.5*(bh_node->ymax - bh_node->ymin);

    if(bh_node->particle != NULL){
        existing_quadrant=get_quadrant(bh_node->particle->x, bh_node->particle->y, bh_node->xmin, bh_node->xmax, bh_node->ymin, bh_node->ymax);

        switch (existing_quadrant) {
            case NE:
                bh_node->NE = create_node(bh_node->particle, xmid, bh_node->xmax, ymid, bh_node->ymax);
                break;
            case NW:
                bh_node->NW = create_node(bh_node->particle, bh_node->xmin, xmid, ymid, bh_node->ymax);
                break;
            case SW:
                bh_node->SW = create_node(bh_node->particle, bh_node->xmin, xmid, bh_node->ymin, ymid);
                break;
            case SE:
                bh_node->SE = create_node(bh_node->particle, xmid, bh_node->xmax, bh_node->ymin, ymid);
                break;
        }

        bh_node->particle = NULL;
    }


    new_quadrant = get_quadrant(new_particle->x, new_particle->y, bh_node->xmin, bh_node->xmax, bh_node->ymin, bh_node->ymax);
    update_mass_of_center(bh_node,new_particle);

    switch (new_quadrant){
        case NE:
            if(bh_node->NE == NULL)
            {
                bh_node->NE = create_node(new_particle, xmid, bh_node->xmax, ymid, bh_node->ymax);
            } else {
                insert_particle(new_particle,bh_node->NE);
            }
            break;
        case NW:
            if(bh_node->NW == NULL)
            {
                bh_node->NW = create_node(new_particle, bh_node->xmin, xmid, ymid, bh_node->ymax);
            } else {
                insert_particle(new_particle,bh_node->NW);
            }
            break;
        case SW:
            if(bh_node->SW == NULL)
            {
                bh_node->SW = create_node(new_particle, bh_node->xmin, xmid, bh_node->ymin, ymid);
            } else {
                insert_particle(new_particle,bh_node->SW);
            }
            break;
        case SE:
            if(bh_node->SE == NULL)
            {
                bh_node->SE = create_node(new_particle, xmid, bh_node->xmax, bh_node->ymin, ymid);
            } else {
                insert_particle(new_particle,bh_node->SE);
            }
            break;
    }
}

// Destroy tree for current delta
void destroy_tree(struct node *bh_node){
  if(bh_node == NULL){
    return;
  }
  if(bh_node->NE != NULL)
      destroy_tree(bh_node->NE);
  if(bh_node->NW != NULL)
      destroy_tree(bh_node->NW);
  if(bh_node->SW != NULL)
      destroy_tree(bh_node->SW);
  if(bh_node->SE != NULL)
      destroy_tree(bh_node->SE);
  free(bh_node);
}

// Update the center mass of the box
void update_mass_of_center(struct node * bh_node, struct part * particle)
{
    bh_node->centerx = update_center(bh_node->centerx, particle->x, bh_node->total_mass, particle->mass);
    bh_node->centerx = update_center(bh_node->centery, particle->y, bh_node->total_mass, particle->mass);
    bh_node->total_mass += particle->mass;
    return;
}

double update_center(double center, double position, double total_mass, double mass){
  return(total_mass * center + mass*position)/(total_mass + mass);
}

// Calculate the force acting on a particle
// Runs the Barnes-Hut algorithm and updates the force
void calculate_forces(struct node *bh_node, struct part *particle, double G, double theta){
    double dx, dy, dist;
    double affecting_mass;
    double force;

    dx = bh_node->centerx - particle->x;
    dy = bh_node->centery - particle->y;
    dist = sqrt(pow(dx,2) + pow(dy,2));

    // check the d/r < theta condition.
    if( (((dist / bh_node->box_size) > theta) || (bh_node->particle))&&(bh_node->particle!=particle) )
    {
        affecting_mass = G * bh_node->total_mass * particle->mass;
        force = affecting_mass / pow(dist,3);;
        particle->fx += force * dx;
        particle->fy += force * dy;
    } else {
        if(bh_node->NE) { calculate_forces(bh_node->NE, particle, G, theta); }
        if(bh_node->NW) { calculate_forces(bh_node->NW, particle, G, theta); }
        if(bh_node->SW) { calculate_forces(bh_node->SW, particle, G, theta); }
        if(bh_node->SE) { calculate_forces(bh_node->SE, particle, G, theta); }
    }
}
