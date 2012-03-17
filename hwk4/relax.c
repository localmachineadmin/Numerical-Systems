#include "mg.h"
#include <math.h>

void relax(float ***u, float ***rhs, int n)
/*
  Red-black Gauss-Seidel relaxation for model problem. Updates the current value of the solution
  u[1..n][1..n], using the right-hand side function rhs[1..n][1..n].
*/
{
int i,j,k;
double dx = 1/(n-1);
double C = alpha * dt/pow(dx,2);
//Jacobi
if (mode==1) 
{
	float ***matnew = u;
	fill0(matnew, n);
	for (i=2; i<n; i++) 
		for (j=2; j<n; j++)
			for (k=2; k<n; k++)	
			{
				matnew[i][j][k] 
				= C/(6*C+1) * (u[i+1][j][k] + u[i-1][j][k] 
				+ u[i][j+1][k] + u[i][j-1][k]
				+ u[i][j][k+1] + u[i][j][k-1])
				- 1/(6*C+1) * rhs[i][j][k];
			}
	copy(u, matnew, n);
} 

//Gauss-Siedel
else if (mode==2)
{
	for (i=2; i<n; i++) 
		for (j=2; j<n; j++)
			for (k=2; k<n; k++)	
			{
				u[i][j][k] 
				= C/(6*C+1) * (u[i+1][j][k] + u[i-1][j][k] 
				+ u[i][j+1][k] + u[i][j-1][k]
				+ u[i][j][k+1] + u[i][j][k-1])
				- 1/(6*C+1) * rhs[i][j][k];
			}
}

//GS with Red-Black Ordering
else if (mode==3)
{
  int i,ipass,isw,j,jsw,k,ksw=1;
  double h,h2;
  h=1.0/(n-1);
  h2=h*h;
  /* Red and black sweeps.*/
  /* jsw and isw toggle between 1 and 2 and
     determine starting row in each column
     for given pass  */
  for (ipass=1;ipass<=2;ipass++,ksw=3-ksw) { 
    jsw=ksw;
	for (k=2;k<n;k++,jsw=3-jsw)
		isw=jsw;
		for (j=2;j<n;j++,isw=3-isw)
      /*Gauss-Seidel formula.*/
			for (i=isw+1;i<n;i+=2)
			{
				u[i][j][k] 
				= C/(6*C+1) * (u[i+1][j][k] + u[i-1][j][k] 
				+ u[i][j+1][k] + u[i][j-1][k]
				+ u[i][j][k+1] + u[i][j][k-1])
				- 1/(6*C+1) * rhs[i][j][k];
			}
  }
}
}
