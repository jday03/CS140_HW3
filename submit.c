/*
Assignment 3 
Team Member 1 :
Team Member 2 :
*/

#include "nBody.h"

void readnbody(double** s, double** v, double* m, int n) {
	int myrank;
	int nprocs;
	int i;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	int size = n / nprocs;
	// Open the input file
	FILE * fp;
	fp = fopen ("./input.txt", "r");
	if (fp == NULL) {
	  fprintf(stderr, "error, cannot find input file");
	}


	tempS = (double **)malloc(sizeof(double *) * size);
	for (i = 0; i < size; i++) {
		tempS[i] = (double*)malloc(sizeof(double) * 3);
		for(j = 0; j < 3; j++) {
			tempS[i][j] = 0;
		}
	}

	TempV = (double **)malloc(sizeof(double *) * size);
	for (i = 0; i < size; i++) {
		TempV[i] = (double*)malloc(sizeof(double) * 3);
		for(j = 0; j < 3; j++) {
			TempV[i][j] = 0;
		}
	}

	TempM = (double *)malloc(sizeof(double) * size);

	for(i = 0; i < size; i++) {
		m[i] = 0;
	}
	
	if (myrank == 0) {
		int i;
		for (i = 0; i < nprocs; i++){
			double x, y, z, vx, vy, vz, mass;

			if (result != 7) {
				fprintf(stderr, "error reading body %d. Check if the number of bodies is correct.\n", i);
				exit(0);
			}
			if( i == 0 ){
				int j;
				for (j = 0; i < size; i++) {
					int result = fscanf(fp, INPUT_BODY, &x, &y, &z, &vx, &vy, &vz, &mass);
					if (result != 7) {
						fprintf(stderr, "error reading body %d. Check if the number of bodies is correct.\n", i);
						exit(0);
					}

					s[j][0] = x;
					s[j][1] = y;
					s[j][2] = z;
					v[j][0] = vx;
					v[j][1] = vy;
					v[j][2] = vz;
					m[j] = mass;
				}
			}
			else{ //MPI_Send to the rest of processors using tempS, tempV, tempM

			}

		}
	}

	else{
		//MPI Receive Data
	}
	fclose(fp);


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

void nbody(double** s, double** v, double* m, int n, int iter, int timestep) {
	int myrank;
	int nprocs;
	int i;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	// This is an example of printing the body parameters to the stderr. Your code should print out the final body parameters
	// in the exact order as the input file. Since we are writing to the stderr in this case, rather than the stdout, make
	// sure you dont add extra debugging statements in stderr.

	if (myrank == 0) {
		for (i = 0; i < n / nprocs; i++) {
			fprintf(stderr, OUTPUT_BODY, s[i][0], s[i][1], s[i][2], v[i][0], v[i][1], v[i][2], m[i]);
		}
	}
}

