#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


/***  ***/


void rWc(double *kap, int *m, double *W)
{
  /*  described in Wood(1994) */
  GetRNGstate();

  double b=( -2.0*(*kap) + sqrt(  4*pow(*kap,2)+pow(*m-1.0,2) ) )/(*m-1.0) ;
  double x0=(1.0-b)/(1.0+b) ;
  double c= (*kap)*x0 +(*m-1.0)*log(1.0-pow(x0,2)) ;
  double Z, U;
  int done=0;

  while( done==0) 
  {
    Z=rbeta( (*m-1.0)/2.0, (*m-1.0)/2.0 );
   *W=( 1-(1+b)*Z)/(1.0-(1.0-b)*Z) ;
    U=runif(0,1);
    if( (*kap)*(*W)+(*m-1)*log(1-x0*(*W))-c  > log(U) )  {done=1;}
  }
  PutRNGstate();
}


/***  ***/


double rtheta_bmf(double k, double a, double b, double c)
{ 
/* target density is of the form  
   dbeta(t,1/2,k)*exp(a*t) * exp(b*sqrt(t))*(exp(c*sqrt(t))+exp(-c*sqrt(t)))

   envelope density is of the form 
   dbeta(t,1/2,g)   where g<=k  */

  double u=2.0;
  double th;
  double g=k; 
  double lrth=0.0;
  double lrmx=0.0;
  double ct;

  if(a>0.0)
  { 
    g = fmax(1.0/(1.0+log(2.0+a)),k-a)     ;
    lrmx = a-k+g + (k-g)*log((k-g)/a) ;
  }
  if(b<=0.0){ lrmx=lrmx+c+log(1.0+exp(-2.0*c))  ; }
  if(b>0)
  { 
    ct=c+log(.5*(1.0+exp(-2.0*c))) ;
    lrmx=lrmx+sqrt(pow(b,2)+pow(ct,2)) + log(2.0) ;
  }

  while(log(u)>lrth-lrmx)
  { 
    u=runif(0,1);
    th=rbeta(.5,g) ;
    lrth=a*th+(k-g)*log(1-th)+
         sqrt(1.0-th)*b+sqrt(th)*c+log(1.0+exp(-2.0*sqrt(th)*c))  ;
   }
  return th; 
}


/***  ***/


void ry_bmfc(double *y, double *l, double *d, int *n)
{
  /* described in Hoff(2009) */

  GetRNGstate();

  int i,j;
  double theta;
  double k= .5*(*n-1.0)  ;
  double smyi, omyi, a, b;

  for(i=0; i<*n; i++)
  { 
    omyi=1.0/(1.0-pow(y[i],2)) ;
    smyi=sqrt(omyi);
    a=l[i]+l[i]*pow(y[i],2)*omyi ;
    b= (-1.0)*y[i]*d[i]*smyi ;
    for(j=0 ; j<*n ; j++)
    {
      a=a-l[j]*pow(y[j],2)*omyi ;
      b=b+y[j]*d[j]*smyi;
    }
    theta=rtheta_bmf(k,a,b,fabs(d[i])) ;
    for(j=0;j<*n;j++){ y[j]=y[j]*sqrt(1.0-theta)*smyi ;}
    y[i]= sqrt(theta)*pow(-1.0,rbinom(1,1.0/(1.0+exp(2.0*sqrt(theta)*d[i])) )) ;
  }
  
  PutRNGstate();
}



/*     */

double rtheta_bing(double k, double a)
{
/* target density is of the form  
   dbeta(t,1/2,k)*exp(a*t) 

   envelope density is of the form 
   dbeta(t,1/2,g)   where g<=k  */

  double u=2.0;
  double th;
  double g=k; 
  double lrth=0.0;
  double lrmx=0.0;

  if(a>0.0)
  { 
    g= fmax(1.0/(1.0+log(2.0+a)),k-a)    ; /* k-a <= g <=k , dec in a */
    lrmx= a - k + g + (k-g)*log((k-g)/a) ;
  }
  
  while(log(u)>lrth-lrmx)
  {
    u=runif(0,1);
    th=rbeta(.5,g);
    lrth=a*th+(k-g)*log(1-th);
  }
  return th;

  }

/*    */

void ry_bingc(double *y, double *l, int *n)
{
  /* described in Hoff(2009) */
  /* should be the same as ry_bmf(y,l,rep(0,n),n) */
  GetRNGstate();

  int i,j;
  double theta;
  double k= .5*(*n-1.0)  ;
  double smyi, omyi, a;

  for(i=0; i<*n; i++)
  {
    omyi=1.0/(1.0-pow(y[i],2)) ;
    smyi=sqrt(omyi);
    a=l[i]+l[i]*pow(y[i],2)*omyi ;
    for(j=0 ; j<*n ; j++)
    {
      a=a-l[j]*pow(y[j],2)*omyi ;
    }
    theta=rtheta_bing(k,a) ;
    for(j=0;j<*n;j++){ y[j]=y[j]*sqrt(1.0-theta)*smyi ;}
    y[i]= sqrt(theta)*pow(-1.0,rbinom(1,0.5)) ;
  }

  PutRNGstate();
}

/*    */

static const R_CMethodDef CEntries[] = {
  {"rWc",      (DL_FUNC) &rWc,      3},
  {"ry_bingc", (DL_FUNC) &ry_bingc, 3},
  {"ry_bmfc",  (DL_FUNC) &ry_bmfc,  4},
  {NULL, NULL, 0}
};

void R_init_rstiefel(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}


