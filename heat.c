//Solving 3D heat in FTCS

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "nrutil.h"

int main(void) 
{
	double lx = 1;
	double ly = 1;
	double lz = 1;
	double alpha = 0.001;
	double dt = 0.005;
	int nsteps = 10;
	double source = 0;
	double noise = 0.99;
	
	//set problem size. set ny=nz=nx for testing purposes
	int nx = 500;
	int ny = nx;
	int nz = nx;
	
	int i, j, k, n;
	time_t start, end;
	double time_diff;
	
	double dx = lx/nx;
	double dy = ly/ny;
	double dz = lz/nz;
	double Cx = alpha*dt/pow(dx,2);

	float ***matT;
	float ***matTnew;
	double *vecX;
	double *vecY;
	double *vecZ;
	
	matT = f3tensor(1, nx, 1, ny, 1, nz);
	matTnew = f3tensor(1, nx, 1, ny, 1, nz);	
	vecX = dvector(1, nx);
	vecY = dvector(1, ny);
	vecZ = dvector(1, nz);
	
	for (i=1; i<=nx; i++)
		vecX[i] = i/nx;
	for (i=1; i<=ny; i++)
		vecY[i] = i/ny;
	for (i=1; i<=nz; i++)
		vecZ[i] = i/nz;		
		
	//establish initial conditions for matrix and random noise
	for (i=1; i<=nx; i++) 
	   for (j=1; j<=ny; j++) 
	      for (k=1; k<=nz; k++) 
		  {
			matT[i][j][k] = exp(
				exp(-pow(5*vecX[i]-2.5,2))*
				exp(-pow(5*vecY[j]-2.5,2))*
				exp(-pow(5*vecZ[k]-2.5,2))
				) 
				* (noise+(rand()%10)/1000);
	      }

	//set initial boundary conditions
	for (i=1; i<=nx; i++) 
		for (j=1; j<=ny; j++)
			for (k=1; k<=nz; k++)
				if (i==1||i==nx+1||j==1||j==ny+1||k==1||k==nz+1)
					matT[i][j][k] = 0;
			
	//do FTCS
	start = clock();
	for (n=1; n<=nsteps; n++)
		for (i=2; i<nx; i++) 
			for (j=2; j<ny; j++)
				for (k=2; k<nz; k++)
				{
					matTnew[i][j][k] = 	matT[i][j][k]
						+ Cx*
						(
						matT[i-1][j][k]+ matT[i+1][j][k]
						+ matT[i][j-1][k]+ matT[i][j+1][k]
						+ matT[i][j][k-1]+ matT[i][j][k+1]
						- 6*matT[i][j][k]
						)
						+ source;
				}
						
	//set boundary
	for (i=1; i<=nx; i++) 
		for (j=1; j<=ny; j++)
			for (k=1; k<=nz; k++)
				if (i==1||i==nx+1||j==1||j==ny+1||k==1||k==nz+1)
					matT[i][j][k] = 0;		 //Currently set as Dirichlet (dt*n) for periodic

	end = clock();
	time_diff = (double)(end-start)/CLOCKS_PER_SEC;
	printf("Problem size: %d, %lf seconds\n",nx, time_diff);
	
	return 0;
}
	
	
		