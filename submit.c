/*
Assignment 3 
Team Member 1 : Noel Vargas
Team Member 2 : Jonathan Day
*/

#include "nBody.h"
#include <time.h>
#define M_PI 3.14159265358979323846

void readnbody(double** s, double** v, double* m, int n) {
    MPI_Status status;
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

    double **tempS;
    int j;
	tempS = (double **)malloc(sizeof(double *) * size);
	for (i = 0; i < size; i++) {
		tempS[i] = (double*)malloc(sizeof(double) * 3);
		for(j = 0; j < 3; j++) {
			tempS[i][j] = 0;
		}
	}
    double ** tempV;
	tempV = (double **)malloc(sizeof(double *) * size);
	for (i = 0; i < size; i++) {
		tempV[i] = (double*)malloc(sizeof(double) * 3);
		for(j = 0; j < 3; j++) {
			tempV[i][j] = 0;
		}
	}

    double * tempM;
	tempM = (double *)malloc(sizeof(double) * size);

	for(i = 0; i < size; i++) {
		m[i] = 0;
	}


	//save to self and send to rest of processors.
	if (myrank == 0) {
		int i;
		for (i = 0; i < nprocs; i++){
			double x, y, z, vx, vy, vz, mass;
    		int result = 7;
			if( i == 0 ){
				int j;
				for (j = 0; j < size; j++) {
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

				for (j = 0; j < size; j++) {
					int result = fscanf(fp, INPUT_BODY, &x, &y, &z, &vx, &vy, &vz, &mass);
					if (result != 7) {
						fprintf(stderr, "error reading body %d. Check if the number of bodies is correct.\n", i);
						exit(0);
					}

					tempS[j][0] = x;
					tempS[j][1] = y;
					tempS[j][2] = z;
					tempV[j][0] = vx;
					tempV[j][1] = vy;
					tempV[j][2] = vz;
					tempM[j] = mass;

				}
				sendItem(tempS, size, 3, i );
				sendItem(tempV, size, 3, i );
                MPI_Send(&tempM[0], size, MPI_DOUBLE, i, 0,MPI_COMM_WORLD );

			}
		}
	}

	else{
		receiveItem(s,size,3,0);
		receiveItem(v,size,3,0);
        MPI_Recv(&m[0],size,MPI_DOUBLE, 0, 0,MPI_COMM_WORLD,&status);

	}

	fclose(fp);
	free(tempM);
	free(tempS);
	free(tempV);

}

void gennbody(double** s, double** v, double* m, int n){
	int myrank;
	int nprocs;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	size = n/nprocs;
	double massMax = 10^30;
	double distMax = 0.5 * (10 ^13);
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
	int i,j;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	//printf("HAVE REACHED NBODY FUNCTION THERE ARE %u", nprocs);

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

	double * mCurrent;
	mCurrent = (double *)malloc(sizeof(double) * sizeOfS);

	for(i = 0; i < sizeOfS; i++) {
		mCurrent[i] = m[i];
	}


	int iterCount;
	for (iterCount = 0; iterCount < iter; ++iterCount) {


        bodiesForceCalc(s, s, f, m, m, sizeOfS, myrank,myrank); // initial calculation within each processor
		int procCount;



if(nprocs != 1){


    int receiveSourceOriginal = myrank - 1;
    if (receiveSourceOriginal < 0) {
        receiveSourceOriginal = nprocs + receiveSourceOriginal;
    }

    int receiveSource = myrank - 1;
    if (receiveSource < 0) {
        receiveSource = nprocs + receiveSource;
    }


    int sendSource = myrank + 1;
    if (sendSource == nprocs) {
        sendSource = 0 + sendSource - nprocs;
    }

        rotateItems(currentS, f,mCurrent, myrank, sizeOfS, sendSource, receiveSource);

        for(procCount = 0; procCount < nprocs - 1; ++procCount) {



            bodiesForceCalc(s, currentS, f, m,mCurrent, sizeOfS, myrank,receiveSourceOriginal);

            receiveSourceOriginal = receiveSourceOriginal - 1;
            if (receiveSourceOriginal < 0) {
                receiveSourceOriginal = nprocs + receiveSourceOriginal;
            }


			rotateItems(currentS, f,mCurrent, myrank, sizeOfS, sendSource, receiveSource);


        }
}
			int tally;
			for (tally = 0; tally < sizeOfS; ++tally) {

                if(isfinite(f[tally][0]) && isfinite(f[tally][1] ) && isfinite(f[tally][2])) {

                    v[tally][0] += (-(f[tally][0] )/ m[tally]) * timestep;
                    v[tally][1] += (-(f[tally][1] )/ m[tally]) * timestep;
                    v[tally][2] += (-(f[tally][2] )/ m[tally]) * timestep;

                    s[tally][0] += v[tally][0] * timestep;
                    s[tally][1] += v[tally][1] * timestep;
                    s[tally][2] += v[tally][2] * timestep;
                }
			}


        zeroFValues(f,sizeOfS);

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



void zeroFValues(double ** f, int sizeOf){
    int iCount;
    int jCount;
    for(iCount = 0; iCount < sizeOf; ++iCount){
        for(jCount=0;jCount < 3;++jCount){
            f[iCount][jCount] = 0;
        }

    }
}



void rotateItems(double ** s, double ** f,double * mCurrent,int myrank, int sizeOfS, int sendSource, int receiveSource){
    int i,j;
	MPI_Status status;

	if (myrank % 2 == 0) {
        sendItem(s, sizeOfS, 3, sendSource);
        sendItem(f, sizeOfS, 3, sendSource);
		MPI_Send(&mCurrent[0],sizeOfS,MPI_DOUBLE, sendSource, 0,MPI_COMM_WORLD );

		receiveItem(s, sizeOfS, 3, receiveSource);
        receiveItem(f, sizeOfS, 3, receiveSource);
		MPI_Recv(&mCurrent[0],sizeOfS,MPI_DOUBLE, receiveSource, 0,MPI_COMM_WORLD,&status);
    } else {
        //because we are receiving first here we must make buffers


		double **coordBuffer;
		coordBuffer = (double **) malloc(sizeof(double *) * sizeOfS);
		for (i = 0; i < sizeOfS; i++) {
			coordBuffer[i] = (double *) malloc(sizeof(double) * 3);
			for (j = 0; j < 3; j++) {
				coordBuffer[i][j] = s[i][j];
			}

		}
		receiveItem(s, sizeOfS, 3, receiveSource);

		double **forceBuffer;
		forceBuffer = (double **) malloc(sizeof(double *) * sizeOfS);
		for (i = 0; i < sizeOfS; i++) {
			forceBuffer[i] = (double *) malloc(sizeof(double) * 3);
			for (j = 0; j < 3; j++) {
				forceBuffer[i][j] = f[i][j];
			}

		}



		receiveItem(f, sizeOfS, 3, receiveSource);

		double * mBuffer;
		mBuffer = (double *)malloc(sizeof(double) * sizeOfS);

		for(i = 0; i < sizeOfS; i++) {
			mBuffer[i] = mCurrent[i];
		}

		MPI_Recv(&mCurrent[0],sizeOfS,MPI_DOUBLE, receiveSource, 0,MPI_COMM_WORLD,&status);
      //  printf("MASS 2 IS %f", mCurrent[0]);
		sendItem(coordBuffer, sizeOfS, 3, sendSource);
        sendItem(forceBuffer, sizeOfS, 3, sendSource);
		MPI_Send(&mCurrent[0],sizeOfS,MPI_DOUBLE, sendSource, 0,MPI_COMM_WORLD );

        free(coordBuffer);
        free(forceBuffer);
		free(mBuffer);
    }


}



void sendItem(double ** item, int sizeOfFirstPart,int sizeOfSecondPart, int destination){
	//printf("HAVE REACHED send item FUNCTION");

	int count;
	for(count = 0; count < sizeOfFirstPart;++count){
		MPI_Send(&item[count][0],sizeOfSecondPart,MPI_DOUBLE, destination, 0,MPI_COMM_WORLD );
	}

}


void receiveItem(double ** item, int sizeOfFirstPart,int sizeOfSecondPart, int source){
	int count;
    MPI_Status status;
    for(count = 0; count < sizeOfFirstPart;++count){
		MPI_Recv(&item[count][0],sizeOfSecondPart,MPI_DOUBLE, source, 0,MPI_COMM_WORLD,&status );
	}

}

void bodiesForceCalc(double** s1,double ** s2, double** f, double* m,double * mCurrent, int sizeInAProc,int procNumber1,int procNumber2){
	int nCount,otherCount;

	for(nCount = 0; nCount < sizeInAProc; ++nCount) {
		for (otherCount = 0; otherCount< sizeInAProc; ++otherCount) {

				calculateOnItemTwo(s1[nCount][0], s1[nCount][1], s1[nCount][2], m[nCount], s2[otherCount][0],
								   s2[otherCount][1], s2[otherCount][2], mCurrent[otherCount], f[otherCount]);

		}
	}
}







void calculateOnItemTwo(double i1,double j1,double k1,double m1,double i2,double j2,double k2,double m2, double * f2Packet){
	// finding differences
	double diffI = i1-i2;
	double diffJ = j1-j2;
	double diffK = k1-k2;

	//printf ("mass1 is %f \n", m1);
   // printf ("mass2 is %f \n", m2);
    double radius = sqrt((diffI*diffI) + (diffJ * diffJ) + (diffK*diffK));
    if(radius != 0 && isfinite(radius)) {
    const double G = 0.6674; //G VALUE HERE
     double scalarForce = (G * m1 * m2) / (radius * radius);

    // computing F* (ijk - ijk)/r
    diffI = diffI / radius;
    diffJ = diffJ / radius;
    diffK = diffK / radius;

    f2Packet[0] += (diffI * scalarForce);
    f2Packet[1] += (diffJ * scalarForce);
    f2Packet[2] += (diffK * scalarForce);


    }

}


