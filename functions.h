int readfile(char *filename, int Nparticles);

double distance(double xj, double xi, double yj, double yi);   

int sort(int Nparticles, double *d, int *ID);

int construct_triangles(struct particle *part);

int triangle_counter(struct particle *part, double cx, double cy); 

double composite_simpson(struct rings * ring, int Npoints, double a, double b);

int my_interpolation(int Nrings, double *S, double *R, gsl_spline * spline, gsl_interp_accel * acc);

double spline_root( gsl_spline * spline, 
		    gsl_interp_accel * acc, 
		    double alpha, 
		    double a0, 
		    double a1, 
		    double epsilon, 
		    int nmax);

