#include "nBody.h"
#include <time.h>
#define M_PI 3.14159265358979323846

int main(int argc, char *argv[]) {

	MPI_Init(&argc,&argv);
	int myrank;
	int nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int n, iters, timestep;

	double** s; // position in 3D space for each body
	double** v; // velocity in 3D space for each body
	double* m;  // mass of each body
	int size; // # of bodies stored on each proc.

	// arguments: ./nBody r #n #iter #timestep
	// or         ./nBody g #n #iter #timestep
	if (argc != 5) {
		fprintf(stderr, "invalid arguments\n");
		fprintf(stderr, "Run with: ./assn3 [r,g] ");
		fprintf(stderr, " #n #iters #timestep\n");
  		MPI_Finalize();
		exit(0);
	}
	n = atoi(argv[2]);
	size = n / nprocs;
	iters = atoi(argv[3]);
	timestep = atoi(argv[4]);
	
	if(myrank == 0){
		fprintf(stderr,"n = %d, iters = %d, timestep = %d\n", n, iters, timestep);
	}
 
	int i; 
	int j;

	s = (double **)malloc(sizeof(double *) * size);
	for (i = 0; i < size; i++) {
		s[i] = (double*)malloc(sizeof(double) * 3);
		for(j = 0; j < 3; j++) {
			s[i][j] = 0;
		}
	}

	v = (double **)malloc(sizeof(double *) * size);
	for (i = 0; i < size; i++) {
		v[i] = (double*)malloc(sizeof(double) * 3);
		for(j = 0; j < 3; j++) {
			v[i][j] = 0;
		}
	}

	m = (double *)malloc(sizeof(double) * size);

	for(i = 0; i < size; i++) {
		m[i] = 0;
	}

	if (strcmp(argv[1], "r") == 0) {
		readnbody(s, v, m, n); 
	} else {
		gennbody(s, v, m, n);
	}
	/*if(myrank==0) {
		fprintf(stderr,"starting nbody...\n");
	}*/
	nbody(s, v, m, n, iters, timestep);

	for (i = 0; i < size; i++) {
		free(s[i]);
		free(v[i]);
	}

	free(s);
	free(v);
	free(m);

	MPI_Finalize();
	return 0;
}

void gennbody(double** s, double** v, double* m, int n){
	//Write This Func:
	int myrank;
	int nprocs;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	size = n/nprocs;
	double massMax = 10^30;
	double distMax = 0.5 *10 ^13;
	double dist;
	double thetaMax = 2 * M_PI;
	double theta;
	double zMax = 10 ^11;




	double div;


	int i, j;
	srand(time(NULL));
	for(i = 0; i < size; i++){
		//generate random masses from 0 to 10^30
		div = RAND_MAX / massMax;
		m[i] = rand() /div;
		//generate random position and set velocity = 0
		for(j=0; j < 3; j++){
			v[i][j] = 0;  //velocity

			// setting positions
			div = RAND_MAX/ distMax;
			dist = rand() / div;
			div = RAND_MAX / thetaMax;
			theta = rand() / div;

			if(j == 0){//x
				s[i][j] = dist * cos(theta);
			}
			if(j == 1){//y
				s[i][j] = dist * sin(theta);
			}
			if(j == 2){ //z
				div = RAND_MAX/ zMax;
				s[i][j] = zMax/2 + (rand() /div );
			}
		}

	}

}

void readnbody(double** s, double** v, double* m, int n){
	//Write This Func
}

void nbody(double** s, double** v, double* m, int n, int iter, int timestep) {
	//Write This Func
}