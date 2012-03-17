#include "mg.h"
#include "nrutil.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv){
	int i,j,k;
	int n = 257; //set problem size here (33, 65, 129, 257)
	
	int ncycle = 2;  
	double dx = 1/(n-1);
	double C = alpha * dt/pow(dx,2);
	double time_diff;
	

	mode = 1; //set as 1 for Jacobi, 2 for GS, 3 for GS with Black-Red ordering
	
	if (mode == 1)
		printf("Starting Jacobi:\n");
	else if (mode == 2)
		printf("Starting Gauss-Siedel:\n");
	else if (mode == 3)
		printf("Starting Gauss-Siedel w/ Red-Black:\n");
	
	time_t start, end;
	FILE *outfile;
	
	double *vec;
	float ***f;

	vec = dvector(1, n);
	f = f3tensor(1, n, 1, n, 1, n);
	
	for (i=1; i<=n; i++)
		vec[i] = i/n;

	//establish initial conditions for matrix and random noise
	for (i=1; i<=n; i++) 
	   for (j=1; j<=n; j++) 
	      for (k=1; k<=n; k++) 
		  {
			f[i][j][k] = exp(
				exp(-pow(5*vec[i]-2.5,2))*
				exp(-pow(5*vec[i]-2.5,2))*
				exp(-pow(5*vec[i]-2.5,2))) 
				* (noise+(rand()%10)/1000);
	      }
		  
	start = clock();
	mglin(f,n,ncycle);
	end = clock();
	
	time_diff = (double)(end-start)/CLOCKS_PER_SEC;
	printf("%lf seconds\n", time_diff);

	outfile = fopen("soln.dat", "w");
	fwrite(&f[1][1][1],sizeof(double),n*n*n,outfile);
	fclose(outfile);
	free_f3tensor(f,1,n,1,n,1,n);
}
