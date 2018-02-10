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

	int sizeOfS = n / nprocs;

	// initializing f array in the same way as s was in nBody.c
	double **f;
	f = (double **) malloc(sizeof(double *) * sizeOfS);
	for (i = 0; i < sizeOfS; i++) {
		f[i] = (double *) malloc(sizeof(double) * 3);
		for (j = 0; j < 3; j++) {
			f[i][j] = 0;
		}
	}



	// initializing current position array for calculation against s
	double **currentS;
	currentS = (double **) malloc(sizeof(double *) * sizeOfS);

	for (i = 0; i < sizeOfS; i++) {
		currentS[i] = (double *) malloc(sizeof(double) * 3);
		for (j = 0; j < 3; j++) {
			currentS[i][j] = 0;
		}
	}


	int iterCount;
	for (iterCount = 0; iterCount < iter; ++iterCount) {


		initialBodiesForceCalc(s, f, m, sizeOfS, myrank);




		int procCount;
		for(procCount = 0; procCount < nprocs - 1; ++procCount) {


			int receiveSource = myrank - 1;
			if (receiveSource < 0) {
				receiveSource = nprocs + receiveSource;
			}


			int sendSource = myrank + 1;
			if (sendSource > nprocs) {
				sendSource = 0 + sendSource - nprocs;
			}


			if (myrank % 2 == 0) {
				send(s, sizeOfS, 3, sendSource);
				send(f, sizeOfS, 3, sendSource);
				receive(s, sizeOfS, 3, receiveSource);
				receive(f, sizeOfS, 3, receiveSource);
			} else {
				//because we are receiveing first here we must make buffers
				double **forceBuffer;
				forceBuffer = (double **) malloc(sizeof(double *) * sizeOfS);
				for (i = 0; i < sizeOfS; i++) {
					forceBuffer[i] = (double *) malloc(sizeof(double) * 3);
					for (j = 0; j < 3; j++) {
						forceBuffer[i][j] = f[i][j];
					}

				}


				receive(f, sizeOfS, 3, iterCount - 1);

				double **coordBuffer;
				coordBuffer = (double **) malloc(sizeof(double *) * sizeOfS);
				for (i = 0; i < sizeOfS; i++) {
					coordBuffer[i] = (double *) malloc(sizeof(double) * 3);
					for (j = 0; j < 3; j++) {
						coordBuffer[i][j] = currentS[i][j];
					}

				}

				receive(currentS, sizeOfS, 3, receiveSource);


				send(coordBuffer, sizeOfS, 3, sendSource);
				send(forceBuffer, sizeOfS, 3, sendSource);

				free(coordBuffer);
				free(forceBuffer);
			}


			int receiveSourceOriginal = myrank - procCount;
			if (receiveSourceOriginal < 0) {
				receiveSourceOriginal = nprocs + receiveSourceOriginal;
			}

				bodiesForceCalc(s, currentS, f, m, sizeOfS, myrank,	receiveSourceOriginal); // INCLUDE THAT FUNCTION PARAMETER SHIT WHEN SEND/RECEIVE HAS ITS SHIT TOGETHER





			// SEND FORCES TO WHERE THEY BELONG AS f HERE

			int tally;
			for (tally = 0; tally < sizeOfS; ++tally) {
				v[tally][0] += -f[tally][0] / m[myrank * sizeOfS + tally];
				v[tally][1] += -f[tally][1] / m[myrank * sizeOfS + tally];
				v[tally][2] += -f[tally][2] / m[myrank * sizeOfS + tally];

				s[tally][0] += v[tally][0] * timestep;
				s[tally][1] += v[tally][1] * timestep;
				s[tally][2] += v[tally][2] * timestep;

			}

		}
	}






	// This is an example of printing the body parameters to the stderr. Your code should print out the final body parameters
	// in the exact order as the input file. Since we are writing to the stderr in this case, rather than the stdout, make
	// sure you dont add extra debugging statements in stderr.
	if (myrank == 0) {
		for (i = 0; i < n / nprocs; i++) {
			fprintf(stderr, OUTPUT_BODY, s[i][0], s[i][1], s[i][2], v[i][0], v[i][1], v[i][2], m[i]);
		}
	}
}




void sendItem(double ** item, int sizeOfFirstPart,int sizeOfSecondPart, int destination){
	int count;
	for(count = 0; count < sizeOfFirstPart;++count){
		MPI_SEND(item[count][0],sizeOfSecondPart,MPI_DOUBLE, destination, 0,MPI_COMM_WORLD );
	}

}


void receiveItem(double ** item, int sizeOfFirstPart,int sizeOfSecondPart, int source){
	int count;
	for(count = 0; count < sizeOfFirstPart;++count){
		MPI_SEND(item[count][0],sizeOfSecondPart,MPI_DOUBLE, source, 0,MPI_COMM_WORLD );
	}


}

void bodiesForceCalc(double** s1,double ** s2, double** f, double* m, int sizeInAProc,int procNumber1,int procNumber2){
	int nCount,otherCount;

	for(nCount = 0; nCount < n; ++nCount) {
		for (otherCount = 0; otherCount< n; ++otherCount) {

			if(nCount != otherCount) {
				calculateOnItemTwo(s1[nCount][0], s1[nCount][1], s1[nCount][2], m[procNumber1 * sizeInAProc + nCount], s2[otherCount][0],
								   s2[otherCount][1], s2[otherCount][2], m[procNumber2 * sizeInAProc + otherCount], f[otherCount]);
			}
		}
	}
}







void initialBodiesForceCalc(double** s, double** f, double* m, int sizeInAProc,int procNumber){
	int nCount,otherCount;

	for(nCount = 0; nCount < n; ++nCount) {
		for (otherCount = 0; otherCount< n; ++otherCount) {

			if(nCount != otherCount) {
				calculateOnItemTwo(s[nCount][0], s[nCount][1], s[nCount][2], m[procNumber * sizeInAProc + nCount], s[otherCount][0],
								   s[otherCount][1], s[otherCount][2], m[procNumber * sizeInAProc + otherCount], f[otherCount]);
			}
		}
	}
}


void calculateOnItemTwo(double i1,double j1,double k1,double m1,double i2,double j2,double k2,double m2, double * f2Packet){
	// finding differences
	double diffI = i1-i2;
	double diffJ = j1-j2;
	double diffK = k1-k2;

	double radius = sqrt((diffI*diffI) + (diffJ * diffJ) + (diffK*diffK));

	const double G = 6.674 * pow(10,-11); //G VALUE HERE

	double scalarForce = (G * m1*m2)/ (radius * radius);

	// computing F* (ijk - ijk)/r
	diffI = diffI /radius;
	diffJ = diffJ /radius;
	diffK = diffK /radius;

	f2Packet[0] += (diffI * scalarForce);
	f2Packet[1] += (diffJ * scalarForce);
	f2Packet[2] += (diffK * scalarForce);

}


