#include "stdio.h"
#include "stdlib.h"
#include "stdarg.h"
#include "math.h"
#include "assert.h"
#include "string.h"
#include "mt.h"

/* Description: compute the stat matrix according to super delta
   pairing scheme.  Notice that for many statistics, we can explore
   symmetry to save 1/2 computing time.  */
										 
void compute_super_delta_test_stat(GENE_DATA* pdata, int* L, float** T, const int symmetry, FUNC_STAT_PAIR func_stat_pair, const void* extra)
/* This is the super-delta equivalence of compute_test_stat. (in mt.c) */
/*L is an array which contains 0,1,2 ... for specifying class label*/
/*T is a 2d (ngenes x ngenes) matrix of test stats needs to
  return; symmetry is an integer taking value in: 0, 1, -1,
  representing "no symmetry", "symmetry", and "skew symmetry".*/
/*func_stat_pair is a function pointer with the following
  protocol float func_stat_pair(float *Y1, float *Y2, int* L, int
  n, float na, void* extra).  See mt.h for the typedef.*/
{
  int i,j;
  switch (symmetry){
  case 0:                       /* no symmetry. */
    for(i=0;i<pdata->nrow;i++){
        /* T[i]=(*func_stat_pair)(pdata->d[i],pdata->d[j],L,pdata->ncol,pdata->na,extra); */
      for(j=0;j<pdata->nrow;j++){
        T[i][j]=(*func_stat_pair)(pdata->d[i],pdata->d[j],L,pdata->ncol,pdata->na,extra);
      }
    }
  case 1:                       /* symmetric matrix */
    /* compute diagonal elemets first */
    for(i=1;i<pdata->nrow;i++){
      T[i][i]=(*func_stat_pair)(pdata->d[i],pdata->d[j],L,pdata->ncol,pdata->na,extra);
    }
    /* compute lower triangle elems, then map them to the upper triagle.*/
    for(i=1;i<pdata->nrow;i++){
      for(j=0;j<i;j++){
        T[i][j]=(*func_stat_pair)(pdata->d[i],pdata->d[j],L,pdata->ncol,pdata->na,extra);
        T[j][i]=T[i][j];
      }
    }
  case -1:                      /* skew-symmetric matrix */
    /* diagonal elemets must be zeros. */
    for(i=1;i<pdata->nrow;i++){
      T[i][i]=0;
    }
    /* compute lower triangle elems, then map them to the upper triagle.*/
    for(i=1;i<pdata->nrow;i++){
      for(j=0;j<i;j++){
        T[i][j]=(*func_stat_pair)(pdata->d[i],pdata->d[j],L,pdata->ncol,pdata->na,extra);
        T[j][i]=-T[i][j];
      }
    }
  }
}

float delta_subtraction_ttest(const float *Y1, const float *Y2, const int* L,const int n, const float na,const void *extra) 
{
  /* compute delta from two genes Y1, Y2. Here simple subtraction is
     used. */
  float delta[n];
  int i;
  for(i=1;i<n;i++){
    delta[i]=Y1[i]-Y2[i];
  }
  /* Now compute the t-stats */
  float num,denum,res;
  res=two_sample_tstat_num_denum(delta,L,n,na,&num,&denum,extra);
  if(res==NA_FLOAT) return NA_FLOAT;
  return num;
}

float delta_subtraction_nstat(const float *Y1, const float *Y2, const int* L,const int n, const float na,const void *extra) 
{
  /* compute delta from two genes Y1, Y2. Here simple subtraction is
     used. */
  float delta[n];
  int i;
  for(i=1;i<n;i++){
    delta[i]=Y1[i]-Y2[i];
  }
  /* Now compute the N-stats */
  float num,denum,res;
  res=Nstat_num_denum(delta,L,n,na,&num,&denum,extra);
  if(res==NA_FLOAT) return NA_FLOAT;
  return num;
}
