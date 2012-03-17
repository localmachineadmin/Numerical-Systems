void rstrct(double **uc, double **uf, int nc)
/* 
   Half-weighting restriction. nc is the coarse-grid dimension. The fine-grid solution is input in
   uf[1..2*nc-1][1..2*nc-1], the coarse-grid solution is returned in uc[1..nc][1..nc].
*/
{
  int ic,iif,jc,jf,ncc=2*nc-1;
  /* Interior points.*/
  for (jf=3,jc=2;jc<nc;jc++,jf+=2) { 
    for (iif=3,ic=2;ic<nc;ic++,iif+=2) {
      uc[ic][jc]=0.5*uf[iif][jf]+0.125*(uf[iif+1][jf]+uf[iif-1][jf]
					+uf[iif][jf+1]+uf[iif][jf-1]);
    }
  }
  /* Boundary points. */
  for (jc=1,ic=1;ic<=nc;ic++,jc+=2) { 
    uc[ic][1]=uf[jc][1];
    uc[ic][nc]=uf[jc][ncc];
  }
  for (jc=1,ic=1;ic<=nc;ic++,jc+=2) {
    uc[1][ic]=uf[1][jc];
    uc[nc][ic]=uf[ncc][jc];
  }
}
