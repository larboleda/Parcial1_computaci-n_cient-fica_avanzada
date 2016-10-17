#include <stdio.h>
#include <stdlib.h>

 /*global variables*/


struct particle

{
  double pos[2];
  int ngb[3]; //keep the neigs and the particle's id itself (vertices of the triangle)
  int triangles;
};

struct particle *part;


struct rings

{
  double width;
  double density;
  double radius[2]; //keep the external and internal radii of rings
  double r;

};

struct rings *ring;

int Nrings;
int Nparticles,g;
int *id;
double *dc;
