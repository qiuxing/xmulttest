#include "stdio.h"
#include "stdlib.h"
#include "stdarg.h"
#include "math.h"
#include "assert.h"
#include "string.h"
#include "mt.h"

  /* One dimensional N-stat wrapper function. It computes the
     N-statistic from two groups of scalars.  So L1, L2 kernels are
     equivalent.  Missing values are represented by na.  if return ==
     NA_FLOAT, then it has some problems to calculate the N-stat, such
     as the count of one class is less than 2.

     Y: the vector of one gene across experiments.
     n: the number of experiments (slides).
     L: the class labelling of each experiments.
     na: the NA representation of gene values.
     extra: This variable indicates how many biological groups in the study. Since we are not sure about how to use N-stat in an ANOVA setting, the current implementation only deals with two group case. So this variable is not used here.  

     Remarks:
     1. num is the N-stat (squared root version), there is no need to
     devide anything so the denominator is set to be constant 1.0.

     2. There exists an n*log(n) algorithm which depends on qsort().
     Please see ~/my_papers/anova-N-distance/src/N_stat_qsort.c.
     However, this theoretically better algorithm was practically very
     slow.  Only when n>=500 or so it starts outperform the original
     n^2 algorithm.  So I excluded it from the real world application.
   */

float nstat(const float *Y, const int* L,const int n, const float na,const void *extra) 
{
  float num,denum,res;
  res=Nstat_num_denum(Y,L,n,na,&num,&denum,extra);
  if(res==NA_FLOAT) return NA_FLOAT;
  return num;
}

float Nstat_num_denum(const float *Y, const int* L,const int n, const float na,float* num, float* denum,const void* extra) 
{
  /* x,y are vector holders for observations in group 0 and group 1.
     Since dynamic allocation is a royal pain in the ..., I decide to
     alloc full length for both x,y first. */
  float *x, *y;
  x=(float *)Calloc(n,float);
  y=(float *)Calloc(n,float);
  int nx=0, ny=0;               /* length of x, y. */

/* Place holder for sum of distances for x,y groups and the between
   group sum of distances. */
  float WGX=0, WGY=0, BG=0;      

  int i,j,class;
  for (i=0; i<n; i++) {
    if (Y[i]==na)               /* NA are discarded */
       continue;
    class=L[i];
    if(class==0){
      x[nx] = Y[i];
      nx +=1;
    }else{
      y[ny] = Y[i];
      ny +=1;
    };
  }

  for (i=1; i<nx; i++) {
    for (j=0; j<i; j++){
      WGX += fabs(x[i] -x[j]);
    }
  }
  WGX /=(nx*nx*0.5);

  for (i=1; i<ny; i++) {
    for (j=0; j<i; j++){
      WGY += fabs(y[i] -y[j]);
    }
  }
  WGY /=(ny*ny*0.5);

  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      BG += fabs(x[i]-y[j]);
    }
  }
  BG /=(nx*ny*0.5);

  /* eliminate super small numbers which may be numerical errors to
     ensure better numerical stability*/
  if(WGX<EPSILON)
    WGX = 0.0;
  if(WGY<EPSILON)
    WGY = 0.0;
  if(BG<EPSILON)
    BG = 0.0;

  *num = sqrt(BG -WGX -WGY);
  *denum = 1.0;
  Free(x);
  Free(y);
  return 1;
}
