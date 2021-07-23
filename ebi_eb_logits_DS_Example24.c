/* Copyright (C)  2021  Bjarki Eldon 

elambda : compute exact expected branch length spectrum for Examples 2.3 and 2.4 in Durrett and Schweinsberg (2005):
@article{durrett2005coalescent,
  title={A coalescent model for the effect of advantageous mutations on the genealogy of a population},
  author={Durrett, Rick and Schweinsberg, Jason},
  journal={Stochastic processes and their applications},
  volume={115},
  number={10},
  pages={1628--1657},
  year={2005},
  publisher={Elsevier}
}
DOI:10.1016/j.spa.2005.04.009

The algorithm for computing the expected branch length spectrum is described here:
@article{birkner2013statistical,
  title={Statistical properties of the site-frequency spectrum associated with $\Lambda$-coalescents},
  author={Birkner, Matthias and Blath, Jochen and Eldon, Bjarki},
  journal={Genetics},
  volume={195},
  number={3},
  pages={1037--1053},
  year={2013},
  publisher={Oxford University Press}
}


This document and any source code it contains  is distributed under the terms of the GNU General Public Licence (version $\ge 3$).  You
should have received a copy of the licence along with this file (see file COPYING).  


    The source codes  described in this document  are  free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This document and the code it contains   is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.



********************** */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <gsl/gsl_vector.h> 
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_sf_pow_int.h> 
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_sf_elementary.h> 
#include <gsl/gsl_sf_gamma.h> 
#include <gsl/gsl_fit.h> 
#include <gsl/gsl_multifit_nlin.h> 
#include <gsl/gsl_sf_exp.h> 
#include <gsl/gsl_sf_log.h> 
#include <gsl/gsl_sort.h> 
#include <gsl/gsl_statistics_double.h> 
#include <gsl/gsl_integration.h> 
#include <gsl/gsl_errno.h> 

/* example for compiling:
/* $ g++ -Wall -O3 -mtune=corei7 -march=native -DNDEBUG clambdakplusbeta.cpp -lm -lgsl -lgslcblas */
/* for simulator implementing  Kingman + DS coal */
/* print matrix for check */



static void printm( int n, gsl_matrix *X  )
{

  int u,v;
  for( u = 1 ; u <= n ; u++){
    for( v = 1 ; v <= n ; v++){
      printf( "%g ", gsl_matrix_get( X, u, v) ); }
    printf("\n"); }

}


/* return the (incomplete) beta function */ 
static double betafunc(double x, double a, double b )
{
  /* the GSL incomplete beta function is normalised by the complete beta function  */
  return( x < 1. ? gsl_sf_beta_inc( a, b, x) * gsl_sf_beta( a, b) : gsl_sf_beta(a,b)  ); 
}



/* compute the total Beta-coal rate of merging blocks */
static double  lbetank( int n, int j,   const  double a, const double K )
{
  /* compute beta rate ; a is alpha */
  /* n is current number of blocks */
  /* j is number of blocks to merge */
  assert( n > 1);
  assert( j > 1);
  assert( j <= n);

  
  /* double x = ( K > 0. ?  ( a < 1. ?  K/( K + 1.0 +  (cconst*a*pow(K, 1.0 - a)/(1.0 - a))) : (kingman ? K/(K + 1.0) : K / (K +  (1 + pow(2, 1. -a)/(a-1.) ) ) ) ) : 1. ) ; */
  /* in the case a < 1 we have m_N -> 1 + pow(2, 1-a)/(a-1) where  a = alpha_2 >= 2; */
  /* in case a = alpha_1 < 1 take 0 < K < 1; use K as the cutoff point */
  double x = (a < 1.0 ? K/10.0 : (K > 0.0 ?  K / (K +  (1 + pow(2, 1. -a)/(a-1.) )) : 1.0)); 
  return(  gsl_sf_choose( (unsigned)n, (unsigned)j ) *  betafunc( x, ((double)j) - a,  ((double)n) + a -  ((double)j) ) / betafunc(x,  2. - a, a) ) ;
}


/* the multiple-merger part according to example  2.4 in  DS (2005)   */
static double  ds_example_twofour_totalrate( const int n, const int j,  const double c_parameter )
{
  return ( 2.0 * c_parameter * gsl_sf_choose( n, j ) * gsl_sf_beta( (double)j, (double)(n-j+1) ) );
}


/* the multiple merger part according to Example 2.3 in DS (2005) */
static double ds_example_twothree_totalrate( const int n, const int j,  const double cparam,  const double pparam)
{

  assert( powl( (long double)pparam, (long double)(j-2) ) > 0.0);
  assert( powl( (long double)(1.0 - pparam), (long double)(n-j) ) > 0.0 ); 
  return( cparam * gsl_sf_choose( n, j ) * ((double)(powl( (long double)pparam, (long double)(j-2) ) * powl( (long double)(1.0 - pparam), (long double)(n-j) ) ) ) );
}

/* return the jump rate of jumping from i to j blocks */
/* coalescent is Kingman + Cconst * Beta( x;   2-a,a) */
static double qij( int i, int j,  const double a, const double K, const double Cconst, const bool kingman)
{
  /* compute jump rate of block counting process */
  assert( i > 1);
  assert( j < i);
  assert( j > 0);

  /*  ************  Beta-coalescent rate 
   return(  (Cconst*lbetank( i, i-j+1, a, K)) +  ( kingman > 0 ? 1.0 : 0.0) * (i-j == 1 ? gsl_sf_choose(i,2) : 0.) ); 
   *********** */

  /* ***********  DS example 2.4  */
  /* Cconst is now the c parameter in DS 2.4 */
  /* the coalescent rate is  Kingman +  Cconst * Beta( j, n - j + 1) */
    return ( ds_example_twofour_totalrate( i, i -j+1,  Cconst) +  ((kingman > 0 ? 1.0 : 0.0) * (i-j == 1 ? gsl_sf_choose(i,2) : 0.)) );
  /* DS Example 2.3 *** */
  /*   return ( ds_example_twothree_totalrate( i, i -j+1,  Cconst,  a)  +  ((kingman > 0 ? 1.0 : 0.0) * (i-j == 1 ? gsl_sf_choose(i,2) : 0.)) );  */
}


static void QP( int n, const double a, const double K, const double Cconst,  const bool kingman,  gsl_matrix * Q, gsl_matrix *P  )
{

  /* compute matrices  qij and pij */

  int i, j;
  double s, x;
  for( i = 2; i <= n ; i++){
    assert( i <= n);
    s = 0. ;
    for( j = 1; j < i ; j++){
      x = qij( i, j, a, K, Cconst, kingman); 
      s = s +  x;
      gsl_matrix_set( Q, i, j, x);
      gsl_matrix_set( P, i, j, x); }
    gsl_matrix_set( Q, i, i, s);
    for( j = 1; j < i ; j++){
      assert( j < i);
      gsl_matrix_set( P, i, j,  gsl_matrix_get( P, i, j)/s ); } }
}


static void gmatrix( int n,  gsl_matrix * G , gsl_matrix * Q, gsl_matrix * P)
{

  int i, k, m ;
  double s = 0.0 ;

  /* initialise the diagonal */
  for ( i = 2 ; i <= n ; i++){
    gsl_matrix_set( G, i, i, 1./gsl_matrix_get( Q, i,i) );}

  for( i = 3; i <= n ; i++){
    assert( i <= n );
    for( m = 2; m < i; m++){
	s = 0. ;
	for( k = m; k < i ; k++){
	  assert( i <= n);
	  assert( k <= n);
	  assert( m <= n);
	  s = s +  gsl_matrix_get( P, i, k) * gsl_matrix_get( G, k, m) ; }
      gsl_matrix_set( G, i, m, s); } }
}


/* compute the matrix  p^{(n)}[k,b] */
static void pnb( int n, gsl_matrix * lkb, gsl_matrix * G, gsl_matrix * P  )
{

  /* j is nprime from Prop A1 in paper */
  /* lnb is the matrix  p^{(nprime)}[k,b]; used for each fixed k  */
  gsl_matrix * lnb = gsl_matrix_calloc( n+1, n+1);
  
  int  k, b, j, i ;
  double s = 0.0 ;
  gsl_matrix_set( lkb, n, 1, 1.0 ); 
  for ( k = 2 ; k < n ; k++){
    for( i = k ; i <= n ; i++){
    for ( b = 1; b <= i - k + 1 ; b++){
      gsl_matrix_set( lnb, i, b, ( k == i ? ( b == 1 ? 1.0 : 0.0) : 0.0 ) );
      s = 0. ;
      for( j = k ; j < i ; j++){
	gsl_matrix_set( lnb,  i, b, gsl_matrix_get( lnb, i,b) +  (b > i - j ? (((double)(b - i + j)) * gsl_matrix_get( lnb, j, b - i + j) * (gsl_matrix_get( P, i, j) * gsl_matrix_get( G, j, k) / gsl_matrix_get( G, i, k)) / ((double)j)) : 0.0 ) +  (b < j ? ((((double)(j-b))*gsl_matrix_get( lnb, j, b) *  (gsl_matrix_get( P, i, j) * gsl_matrix_get( G, j, k) / gsl_matrix_get( G, i, k)) / ((double)j) )) : 0.0 ) ) ;
      } } }

    for( j = 1 ; j <= (n - k + 1) ; j++){
      gsl_matrix_set( lkb, k, j, gsl_matrix_get( lnb, n, j) ); }
    gsl_matrix_set_zero( lnb ); }

  gsl_matrix_free( lnb ); 
}


static double ltwoebi( const int n, const  double a, const  double K, const double Cconst, const bool kingman )
{

  gsl_matrix * P = gsl_matrix_calloc( n+1, n+1);
  gsl_matrix * Q = gsl_matrix_calloc( n+1, n+1);
  gsl_matrix * G = gsl_matrix_calloc( n+1, n+1);
  gsl_matrix * Pn = gsl_matrix_calloc( n+1, n+1);

  const double  ice[] = {0, 296,63,28, 8,15, 5, 4, 6, 1, 1, 2, 2, 6, 2, 4, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 2, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 3, 1, 5, 0, 2, 0, 0, 2, 2, 2, 5, 5, 9,16,60};
  const double icesegsites = 574.0; 
  double   omega[] =  {0, 354,67,28,16,10, 2, 5, 2, 2, 4, 3, 3, 5, 6, 4, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 6, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 2, 0, 4, 0, 2, 2, 3, 1, 0, 0, 1, 1, 1, 5, 8,12,24,62};
  double osegsites = 656.0;
    
  /* compute the Q and P matrices */
  QP( n, a, K, Cconst, kingman,   Q, P);
  gmatrix( n,  G, Q, P);
  pnb( n, Pn, G, P);
  
  int b, k;
  double s = 0.0 ;
  double eb = 0.0 ;
  double * ebi  = (double *)calloc( n, sizeof(double) ); 

  for ( b = 1; b < n ; b++){
    s = 0. ;
    for( k = 2; k <= n-b+1 ; k++){
      s = s +  (gsl_matrix_get( Pn, k, b) * ((double)k) * gsl_matrix_get( G, n, k)); }
    ebi[b] = s ; 
    eb = eb + s ;}
  
  double result = 0.0;
  for( b = 1; b < n ; b++){
    /* printf( "%g\n", ebi[b]/eb ); } */
    result +=  pow( (ebi[b]/eb) - (omega[b]/osegsites), 2.0) ; }

  gsl_matrix_free( P);
  gsl_matrix_free( Q);
  gsl_matrix_free( G);
  gsl_matrix_free( Pn);
  free( ebi);

  return( sqrt( result) );
}


static void ebi(  int n, const  double a, const  double K, const double Cconst, const bool kingman )
{  
  gsl_matrix * P = gsl_matrix_calloc( n+1, n+1);
  gsl_matrix * Q = gsl_matrix_calloc( n+1, n+1);
  gsl_matrix * G = gsl_matrix_calloc( n+1, n+1);
  gsl_matrix * Pn = gsl_matrix_calloc( n+1, n+1);
    
  /* compute the Q and P matrices */
  QP( n, a, K, Cconst, kingman,   Q, P);
  gmatrix( n,  G, Q, P);
  pnb( n, Pn, G, P);
  
  int b, k;
  double s = 0.0 ;
  double eb = 0.0 ;
  double * ebi  = (double *)calloc( n, sizeof(double) ); 

  for ( b = 1; b < n ; b++){
    s = 0. ;
    for( k = 2; k <= n-b+1 ; k++){
      s = s +  (gsl_matrix_get( Pn, k, b) * ((double)k) * gsl_matrix_get( G, n, k)); }
    ebi[b] = s ; 
    eb = eb + s ;}

  double hali = 0.0; 
  for ( b = 1; b < n ; b++){
    /*  hali += ebi[b]; } */
    printf( "%g\n", log( ebi[b]/eb )  - log( 1.0 -  ebi[b]/eb ) )  ; }
  /*  printf("%g %g\n", ebi[1]/eb, hali/eb) ; */

  gsl_matrix_free( P);
  gsl_matrix_free( Q);
  gsl_matrix_free( G);
  gsl_matrix_free( Pn);
  free( ebi);
}


static void  matun( )
{

  double l, minl,  maxa, maxk, maxc, maxone;
  l = 0.0; maxa = 0.0; maxk = 0.0; minl = 100000.0; maxc = 0.0; maxone = 0.0 ;
  /* sample size, alpha, K, Cconst */

  const int samplesize = 68; 
  for (double one = 0.0; one <= 0.0; one++){
  for( double cconst = 1; cconst <= 10000; cconst*=10.0){
  for ( double a = 80; a < 99; a++){
    for( double k = 1; k < 10; k++){
      /* static double ebi( int n, const  double a, const  double K, const double Cconst, const double kingman ) */
	l = ltwoebi( samplesize, one + a/100.0,  k/10.0,  cconst,   1 );
      maxa = (l < minl ? one + a/100.0 : maxa );
      maxk = (l < minl ? k/10.0 : maxk );
      maxc = (l < minl ?  cconst : maxc );
      maxone = (l < minl ?  one : maxone );
      minl = (l < minl ? l : minl); }}}}

  printf( "%g %g %g %g %g\n", maxone,   maxa, maxk, maxc,  minl);
  ebi( samplesize, maxa, maxk, maxc, 1); 
  /* 
  int a; 
  for( a = 1; a < 100; a++){  
    ebi( atoi(argv[1]), 1. + (((double)a)/100.0),  atof(argv[3]) ); } */
}

int main( int argc, char *argv[])
{
  /* static void ebi(  int n, const  double a, const  double K, const double Cconst, const bool kingman ) */
  /*  ebi( atoi(argv[1]), atof(argv[2])/100.0,  atof(argv[3]), atof(argv[4]), 1); */
  /* DS Example 2.4 */
  /* sample size; c parameter */
   ebi( atoi(argv[1]), 0.01, 0.0, atof(argv[2]), 1);

  /* DS Example 2.3 -  */
  /* sample size; p parameter;  c parameter */
  /*  ebi( atoi(argv[1]), atof(argv[2]), 0.0, atof(argv[3]), 1); */

    
  return GSL_SUCCESS ; 
}
