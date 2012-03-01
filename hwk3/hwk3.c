//Solving 3D heat in Jacobi, Gauss-Siedel and SOR
//compile with gcc -lm -o -hwk3 hwk3.c nrutil.c nrutil.h
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "nrutil.h"

//can set these to specific values
#define breakpoint 0.000001
#define nsteps 100
#define MAX_ITER 1000
#define alpha 0.001
#define dt 0.005
#define noise 0.99

void Jacobi(float ***matT, double Cx, double Cy, double Cz, int nx, int ny, int nz);
void GSiedel(float ***matT, double Cx, double Cy, double Cz, int nx, int ny, int nz);
void SOR(float ***matT, double Cx, double Cy, double Cz, int nx, int ny, int nz);
void duplicate(float ***mat1, float ***mat2, int nx, int ny, int nz);
void dirichlet(float ***mat, int nx, int ny, int nz);


int main(void) 
{
	int nx;
	int i,j,k;
	//set problem size. set ny=nz=nx for testing purpose:
	printf("set problem size (nx=ny=nz):\n");
	scanf("%d", &nx);
	int ny = nx;
	int nz = nx;
	
	double lx = 1;
	double ly = 1;
	double lz = 1;
	double dx = lx/nx;
	double dy = ly/ny;
	double dz = lz/nz;
	double Cx = alpha*dt/pow(dx,2);
	double Cy = alpha*dt/pow(dy,2);
	double Cz = alpha*dt/pow(dz,2);

	float ***matT;

	double *vecX;
	double *vecY;
	double *vecZ;
	
	matT = f3tensor(1, nx, 1, ny, 1, nz);

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

	//set initial dirichlet boundary conditions
	dirichlet(matT, nx, ny, nz);

	//Main methods: comment the other ones out:
	//Jacobi(matT, Cx, Cy, Cz, nx, ny, nz);
	//GSiedel(matT, Cx, Cy, Cz, nx, ny, nz);
	SOR(matT, Cx, Cy, Cz, nx, ny, nz);
	
	return 0;
}
	
//do Jacobi
void Jacobi(float ***matT, double Cx, double Cy, double Cz, int nx, int ny, int nz)
{
	float ***matTnew = f3tensor(1, nx, 1, ny, 1, nz);
	float ***matTold = f3tensor(1, nx, 1, ny, 1, nz);

	int i, j, k, n, iterN;
	time_t start, end;
	double diff = 0;
	double time_diff;
	dirichlet(matTnew, nx, ny, nz);

	for (n=1; n<=nsteps; n++)
	{
		start = clock();
		duplicate(matT, matTold, nx, ny, nz);
		for (iterN = 1; iterN <= MAX_ITER; iterN++)
		{
			for (i=2; i<nx; i++) 
				for (j=2; j<ny; j++)
					for (k=2; k<nz; k++)	
					{
						matTnew[i][j][k] 
						= Cx/(6*Cx+1) * (matT[i+1][j][k] + matT[i-1][j][k])
						+ Cy/(6*Cy+1) * (matT[i][j+1][k] + matT[i][j-1][k])
						+ Cy/(6*Cy+1) * (matT[i][j][k+1] + matT[i][j][k-1])
						+ 1/(2*Cx+2*Cy+2*Cz+1) * matTold[i][j][k];
					}
			//Check for breakpoint if we've surpasses threshold
			diff = 0;
			for (i=1; i<nx; i++) 
				for (j=1; j<ny; j++)
					for (k=1; k<nz; k++)	
						diff += fabs(matT[i][j][k] - matTnew[i][j][k]);
			if (diff/(nx+ny+nz) < breakpoint)
				break;
			duplicate(matT, matTnew, nx, ny, nz);
		}
		if (iterN==MAX_ITER)
			printf("Not converging after %d iterations", iterN);
		end = clock();
		time_diff = (double)(end-start)/CLOCKS_PER_SEC;
		printf("Jacobi Problem size: %d, timestep: %d, %lf seconds\n",nx, n, time_diff);
	}
}

//Do Gauss-Siedel
void GSiedel(float ***matT, double Cx, double Cy, double Cz, int nx, int ny, int nz)
{
	float ***matTlast = f3tensor(1, nx, 1, ny, 1, nz);
	float ***matTold = f3tensor(1, nx, 1, ny, 1, nz);

	int i, j, k, n, iterN;
	time_t start, end;
	double diff = 0;
	double time_diff;

	for (n=1; n<=nsteps; n++)
	{
		start = clock();
		duplicate(matTold, matT, nx, ny, nz);
		for (iterN = 1; iterN <= MAX_ITER; iterN++)
		{
			duplicate(matTlast, matT, nx, ny, nz);
			for (i=2; i<nx; i++) 
				for (j=2; j<ny; j++)
					for (k=2; k<nz; k++)	
					{
						matT[i][j][k] 
						= Cx/(6*Cx+1) * (matT[i+1][j][k] + matT[i-1][j][k])
						+ Cy/(6*Cy+1) * (matT[i][j+1][k] + matT[i][j-1][k])
						+ Cy/(6*Cy+1) * (matT[i][j][k+1] + matT[i][j][k-1])
						+ 1/(2*Cx+2*Cy+2*Cz+1) * matTold[i][j][k];
					}
			//Check for breakpoint if we've surpasses threshold
			diff = 0;
			for (i=1; i<nx; i++) 
				for (j=1; j<ny; j++)
					for (k=1; k<nz; k++)	
						diff += fabs(matT[i][j][k] - matTlast[i][j][k]);
			if (diff/(nx+ny+nz) < breakpoint)
				break;
		}

		if (iterN==MAX_ITER)
			printf("Not converging after %d iterations", iterN);
		end = clock(); 
		time_diff = (double)(end-start)/CLOCKS_PER_SEC; //display timer for timestep
		printf("GS Problem size: %d, timestep: %d, %lf seconds\n",nx, n, time_diff);
	}
}

//Do Successive Over-relaxation
void SOR(float ***matT, double Cx, double Cy, double Cz, int nx, int ny, int nz)
{
	float ***matTlast = f3tensor(1, nx, 1, ny, 1, nz);
	float ***matTold = f3tensor(1, nx, 1, ny, 1, nz);

	double w = 1.65;
	int i, j, k, n, iterN;
	time_t start, end;
	double diff = 0;
	double time_diff;

	for (n=1; n<=nsteps; n++)
	{
		start = clock();
		duplicate(matTold, matT, nx, ny, nz);
		for (iterN = 1; iterN <= MAX_ITER; iterN++)
		{
			duplicate(matTlast, matT, nx, ny, nz);
			for (i=2; i<nx; i++) 
				for (j=2; j<ny; j++)
					for (k=2; k<nz; k++)	
					{
						matT[i][j][k] 
						= ((1-w)*matT[i][j][k]) //This particular line (in particular the (1-w)*matT snippet
												//vastly increases the overall time for SOR to run,
												//by a magnitude of about 100 compared to Jacobi
												//or Gauss-siedel. I am unsure of how to fix this
						+ w*Cx/(6*Cx+1) * (matT[i+1][j][k] + matT[i-1][j][k])
						+ w*Cy/(6*Cy+1) * (matT[i][j+1][k] + matT[i][j-1][k])
						+ w*Cy/(6*Cy+1) * (matT[i][j][k+1] + matT[i][j][k-1])
						+ w/(2*Cx+2*Cy+2*Cz+1) * matTold[i][j][k];

					}
			//Check for breakpoint if we've surpasses threshold
			diff = 0;
			for (i=1; i<nx; i++) 
				for (j=1; j<ny; j++)
					for (k=1; k<nz; k++)	
						diff += fabs(matT[i][j][k] - matTlast[i][j][k]);
			if (diff/(nx+ny+nz) < breakpoint)
				break;
		}

		if (iterN==MAX_ITER)
			printf("Not converging after %d iterations", iterN);
		end = clock(); 
		time_diff = (double)(end-start)/CLOCKS_PER_SEC; //display timer for timestep
		printf("SOR Problem size: %d, timestep: %d, %lf seconds\n",nx, n, time_diff);

	}
}

//copying one matrix into another
void duplicate(float ***matdest, float ***matsource, int nx, int ny, int nz)
{
	int i, j, k;
	for (i=1; i<nx; i++) 
		for (j=1; j<ny; j++)
			for (k=1; k<nz; k++)
				matdest[i][j][k] = matsource[i][j][k];
}

//dirichlet boundary maker
void dirichlet(float ***mat, int nx, int ny, int nz)
{
	int i, j, k;
	for (i=1; i<=nx; i++) 
		for (j=1; j<=ny; j++)
			for (k=1; k<=nz; k++)
				if (i==1||i==nx+1||j==1||j==ny+1||k==1||k==nz+1)
					mat[i][j][k] = 0;
}



	
		