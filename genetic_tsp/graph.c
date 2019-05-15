#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "structs.h"

// populate population into graph
int load_graph() {
  int i, id;
  double x, y;
  char * line;
  char buffer[100];
  int number_of_cities;
  FILE *input_file = fopen(FILENAME, "r");
  line = fgets(buffer, 100, input_file); // NAME
  line = fgets(buffer, 100, input_file); // TYPE
  line = fgets(buffer, 100, input_file); // COMMENT
  line = fgets(buffer, 100, input_file); // DIMENSION
  sscanf(line,"DIMENSION : %d", &number_of_cities);
  line = fgets(buffer, 100, input_file); // EDGE_WEIGHT_TYPE
  line = fgets(buffer, 100, input_file); // NODE_COORD_SECTION
  graph.number_of_cities = number_of_cities;

  for(i=0; i < number_of_cities; i++) {
    line = fgets(buffer, 100, input_file);
    sscanf(line, "%d %lf %lf", &id, &x, &y);
    graph.nodes[i].id = id;
    graph.nodes[i].x = x;
    graph.nodes[i].y = y;
  }

  fclose(input_file);
  return 0;
}
