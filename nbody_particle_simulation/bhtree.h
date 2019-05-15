struct node {
    double total_mass;
    double centerx, centery;
    double xmin, xmax;
    double ymin, ymax;
    double box_size;
    struct part * particle;
    struct node * NE;
    struct node * NW;
    struct node * SW;
    struct node * SE;
};
