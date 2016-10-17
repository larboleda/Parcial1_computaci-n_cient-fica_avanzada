#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "allvars.c"

/*routine to read the file with the particle coordinates.*/

int readfile(char *filename, int Nparticles)
{
  int i;
  FILE *pf = NULL;
  pf = fopen(filename,"r");

  for(i=0;i<Nparticles;i++)

    {
      fscanf(pf,"%lf %lf\n",&part[i].pos[0],&part[i].pos[1]);
      
    }
  
  fclose(pf);

  return 0;
}

/*distance between two particles*/

double distance(double xj, double xi, double yj, double yi)
{

  double d;
  d = sqrt((xj-xi)*(xj-xi) +(yj-yi)*(yj-yi));
  
  return d;
}

/*sorting indexed arrays..*/

int sort(int Nparticles, double *d, int *ID)
{
  
  int i,j, tempID;
  float temp;
  int N;
  N = Nparticles -1;

  for(i=0;i<(N-1);i++)
    
    {
      for(j=0;j< N -1 -i; ++j)
	{
	  if(d[j] > d[j+1])
	    {
	      
	      temp = d[j+1];
	      d[j+1] = d[j];
	      d[j] = temp;
	      
	      tempID = ID[j+i];
	      ID[j+1] = ID[j];
	      ID[j] = tempID;
	      
	    }
	}	      
    }

  return 0;
}


/*routine to look the two nearest neighbors of each particle. It need the positions of the partciles */             

int construct_triangles(struct particle *part)
{
  FILE *pf = NULL;  
  int i,k,j;

  double *d = NULL;
  size_t *ID = NULL;

  d = (double *)malloc((size_t)(Nparticles-1) * sizeof(double));                                                                                                    
  ID = (size_t *)malloc((size_t)(Nparticles-1) *sizeof(size_t));                                                                                                     
  
  pf = fopen("vertices.dat","w");                                                               
                                                                                                                                                                     
  for(i=0;i<Nparticles;i++)                                                                                                                                            
    {                                                                                                                                                                  
      k = 0;                                                                                                                                                           
      id[i]=i;                                                                                                                                               
                                                                                                                                                                      
      for(j=0;j<Nparticles;j++)                                                                                                                                        
        {                                                                                                                                                              
          if(i != j)                                                                                                                                                   
            {                                                                                                                                                          
              d[k]=distance(part[j].pos[0],part[i].pos[0],part[j].pos[1],part[i].pos[1]);                                                                             
              k++;                                                                                                                                                     
            }
	}
      //sorting the array: if you want to use bubble sort algorithm, call the funtion "sort" as: sort(Nparticles,d,ID);                                          
      //this time we will use gsl to made the sorting        
      gsl_sort_index(ID,d,1,Nparticles-1);                                                                                                                             
                                                                                                                                                                       
      //add neighbors of the particle i to the structure                                                                                                               
      part[i].ngb[0] = id[i];                                                                                                                                          
      part[i].ngb[1] = ID[0];  
      part[i].ngb[2] = ID[1];                                                                                                                                          
                                                                                                                                                                   
      fprintf(pf,"%d %d %d\n", part[i].ngb[0], part[i].ngb[1], part[i].ngb[2]); //vertices of the triangle                                                             
    }

  printf("the file with the vertices was wrtten\n");
  
  free(ID);      
                                                                                                                                                                 
  fclose(pf);                                                                                                                                                          

  return 0;
}

/*routine to count the number of triangles formed by a single particle
It needs to compute the from every particle to the center, so cx and cy are the center coordinates.
*/
int triangle_counter(struct particle *part, double cx, double cy) 
{                                                           

  int *v1 = NULL, *v2 = NULL;
  FILE *pf = NULL;
  int i,j;

  v1 = (int *)malloc((size_t)Nparticles * sizeof(int));                                                                                                               
  v2 = (int *)malloc((size_t)Nparticles * sizeof(int));                            
 
  pf = fopen("vertices.dat","r");                                                                                                                                      

  printf("counting triangles..\n");                                                                                                                                    

  for(j=0;j<Nparticles;j++)                                                                                                                                            
                                                                                                                                                                       
    {                                                                                                                                                                  
      //distance from  particles to the center                                                                                                                         
      
      dc[j]=distance(part[j].pos[0], cx, part[j].pos[1], cy);                                                                                                          
      part[j].triangles = 1;                                                                                                                                           
      
      for(i=0;i<Nparticles;i++)                                                                                                                                       
        {                                                                                                                                                              
          fscanf(pf,"%d %d %d\n",&(id[i]),&(v1[i]),&(v2[i]));                                                                                                          
                                                                                                                                                                       
          if(v1[i] == j)                                                                                                                                               
            {
              part[j].triangles += 1;
            }                           

	  if(v2[i] == j)                                                                                                                                               
            {                                                                                                                                                          
             part[j].triangles += 1 ;                                                                                                                                 
           }                                                                                                                                                          
	}                                                                                                                                                              
      
    }                                                                                                                                                  
  free(v1);                                                                                                                                                            
  free(v2);          
  
  fclose(pf);

//write a file with the data of distances to the center and number of triangles per particle
  pf = fopen("num_triang_distances.dat","w");
  printf("written file with distances and number of triangles..\n");
  
  for(i=0;i<Nparticles;i++)                                                                                                                                            
    {                                                                                                                                                                  
      fprintf(pf,"%lf %d\n",dc[i],part[i].triangles);                                                                                                                  
    }                                                                                                                                                        
  
  fclose(pf);

  gsl_sort(dc,1,Nparticles); //organice the dc array in ascending numerical order to assing particles in the rings 


  return 0;
}

/* composite simpson rule to computee the total mass inside a radii r, computed as the integral of the surface density over the ring's area
 a and b are the end points of the interval.
"struct rings *ring" contains the information of the density and the radii of each ring
*/

double composite_simpson(struct rings * ring, int Npoints, double a, double b)  
{

  int i;
  double h,c,P,I;
  double Mtot;

  c = 0.33333333333;
  h = (b-a)/ Npoints; //a and b are the interval ends, in this case b = Rf

  P = 0.0;
  I = 0.0;  

  if (Npoints == 0)
    Mtot = 0;
  else
    {
      for(i=0;i<=Npoints/2 -1;i++)
	
	{
	  P += ring[2*i].density * ring[2*i].r;
	  I += ring[2*i +1].density * ring[2*i + 1].r;
	}
            
      Mtot = 2 * M_PI * h * c *(4.0*I + 2.0*P + ring[Npoints-1].density * b);
    }

  return Mtot;
      
}


/*routine to interpolate the points: surface density as a function of r*/

int my_interpolation(int Nrings, double *S, double *R, gsl_spline * spline, gsl_interp_accel *acc)

{
  
  double Ri, Si;
  FILE *file = NULL;
  file = fopen("interpolation.dat","w");
  
  //gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, (size_t) Nrings);

  gsl_spline_init (spline, R, S, (size_t) Nrings);
 //gsl_interp_accel *acc = gsl_interp_accel_alloc (); //returns a pointer to an accelerator object                                                                    
                 
  for (Ri = R[0]; Ri < R[Nrings -1]; Ri += 0.01)
    {
      Si = gsl_spline_eval (spline, Ri, acc);
      fprintf (file,"%lf\t %lf\n", Ri, Si);
    }
  
  //gsl_spline_free (spline);
  //gsl_interp_accel_free (acc);
  
  
  fclose(file);
  
  return 0;
}


/*routine to find the roots. It returns Re at which:  density = (central_density / e) */

double spline_root( gsl_spline * spline, gsl_interp_accel * acc, double alpha, double a0, double a1, double epsilon, int nmax)
{

  double a2, err, fa0, fa1, fa2;
  int nsteps = 0;

  fa0 = gsl_spline_eval(spline,a0,acc) - alpha;
  fa1 = gsl_spline_eval(spline,a1,acc) - alpha;

  printf("finding root of %e\n",alpha);
  printf("a0 = %lf \t a1 = %lf\n",a0,a1);
  printf("\t fa0 = %lf \t fa1 = %lf\n",fa0,fa1);


  if(a1 > a0)
  
    {
  
      if (fa0*fa1 > 0)
	{
	  printf("warning: function has same sgn on both bounds\n");
	  printf("root may not be found or may not be unique\n");
	  printf("proceed with caution\n");
	}
      
      a2 = 0.5*(a1+a0);
      err = 0.5*fabs(a0-a1);
  
      fa2 = gsl_spline_eval(spline,a2,acc) - alpha;

      
      while (err > epsilon && nsteps < nmax){

	if (fa0*fa2 < 0){
      
	  //a0 = a0
	  a1 = a2;
	  fa1 = gsl_spline_eval(spline,a1,acc) - alpha;
	  
	}
	
	else{
	  
	  //a1 = a1
	  a0 = a2;
	  fa0 =  gsl_spline_eval(spline,a0,acc) - alpha;
	  
	}
	
	++nsteps;
	a2 = 0.5*(a1+a0);
	err = 0.5*fabs(a0-a1);
	
	fa2 = gsl_spline_eval(spline,a2,acc)-alpha; //we must shift the fuction by alpha in order to use the bisection method
      }
      
      if (nsteps >= nmax){
	printf("Did not converge. Do not trust this answer.\n");
      }
 
    }
     return a2;
}


