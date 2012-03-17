#include "mg.h"

void slvsml(float ***u, float ***rhs)
/* 
   Solution of the model problem on the coarsest grid, where h = 1
   2 . The right-hand side is input
   in rhs[1..3][1..3] and the solution is returned in u[1..3][1..3].
*/
{
  void fill0(float ***u, int n);
  double h=0.5;
  double C = alpha/(h*h);
  fill0(u,3);
  u[2][2][2] = -rhs[2][2][2]/(6.0*C+1);
}

//done