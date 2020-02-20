/*
!*****************************************************************************
!                                                                            *
!  Module 'hratelib' for computing heating rates in expanding ejecta         *
!  Original version : OK 29.08.2019                                          *
!  Translated to c++ : HL Dec.23.2019                                        *
!                                                                            *
!*****************************************************************************
*/

#ifndef hratelib_h
#define hratelib_h

#include <vector>
#include <iostream>
#include <cmath>
#include <array>

#include "params.h"
#include "body.h"
#include "tree.h"
#include "utils.h"
#include "space_vector.h"

#define M_PI 3.14159265358979323846  /* pi */

// grid of velocity and Ye from which approximant is interpolated
// All values are for ejecta mass = 0.01 M_/odot
   double YE_GRID[]= {.05, .10, .15, .20, .25, .3, .35, .40, .45, .5};
   double  V_GRID[]= {.05, .1, .2, .3, .4, .5};

// approximant coefficients on the grid
   
   int V_GRID_LEN = std::size(V_GRID);
   int YE_GRID_LEN = std::size(YE_GRID);

   double E0_GRID[YE_GRID_LEN][V_GRID_LEN] = {{5.81, 6.76,  6.5,  6.6,  6.6,  6.6},
                                              { 6.6,  6.6,  9.8,  9.8,  9.8,  9.8}, 
                                              { 6.5,   9.,  10.,  10.,  10.,  10.}, 
                                              { 6.5,  6.5,  5.8, 19.2,  18.,  18.}, 
                                              { 8.4, 29.7, 45.2, 33.2, 33.2,  48.}, 
                                              { 6.1,  18., 47.1, 47.1, 74.8, 74.8}, 
                                              { 7.3,  7.3, 16.3, 23.2, 43.2, 150.}, 
                                              {  .1,  .08,  .03,   .1,   .5,  .15}, 
                                              { 1.2,  1.4, 1.68, 1.68, 1.68, 1.68}, 
                                              {  .4, 1.07,  .76, 1.42, 1.42, 1.42}};

   double ALP_GRID[YE_GRID_LEN][V_GRID_LEN] = {{ 1.35,  1.35,  1.35,  1.35,  1.35,  1.35}, 
                                               { 1.35,  1.35,  1.35,  1.35,  1.35,  1.35}, 
                                               { 1.32,  1.34,  1.36,  1.32,  1.35,  1.35}, 
                                               { 1.36,  1.36,  1.34,  1.34,  1.34,  1.34}, 
                                               {1.368,  1.34, 1.394,  1.8,  1.66,   1.68}, 
                                               { 1.33,  1.33,  1.33,  1.33, 1.374, 1.374},
                                               {1.268, 1.358, 1.384, 1.384, 1.384, 1.344},
                                               { 1.28,  1.28,  1.23,  1.23,   1.1,  1.05}, 
                                               {   8.,    8.,   5.8,   5.8,   5.8,   5.8},  
                                               {  1.1,   1.1,   1.1,   1.1,   1.1,   1.1}};

   double T0_GRID[YE_GRID_LEN][V_GRID_LEN] = {{1.695, 1.34, 1.20, 1.20, 1.20, 1.20}, 
                                              {1.395,1.080, .825, .825, .785, .785}, 
                                              {1.055,  .74,  .64,  .61,  .61,  .61}, 
                                              {  .83,  .60, .455,   .4,  .38,  .38}, 
                                              { .645, .355,  .22,  .17,  .14,  .10}, 
                                              { .540,  .29,  .18,  .13, .095, .081}, 
                                              { .385, .235,   .1,  .06, .035, .025}, 
                                              {  .01,  .01,    0, .023, .013, .013}, 
                                              { .195, .105,  .05,  .03, .023, .021}, 
                                              { .128, .058,  .05, .024,  .02,  .02}};

   double SIG_GRID[YE_GRID_LEN][V_GRID_LEN] = {{  .1,   .1,   .1,   .1,   .1,   .1}, 
                                               {  .1,   .1,  .08,  .08, .075, .075}, 
                                               {.075, .075,  .09,  .06,  .08,  .08}, 
                                               {.075, .075,  .07, .035,  .04,  .04}, 
                                               { .07, .021,  .02, .067, .074, .074}, 
                                               {.096, .047, .021, .021, .017, .017}, 
                                               {.058, .094, .068,  .05,  .03,  .01}, 
                                               {.132, .132,  .05, .004,.0015,.0015}, 
                                               {.144, .091, .026, .019, .015, .012}, 
                                               {.018, .004, .004,.0024,.0018,.0014}};

   double ALP1_GRID[YE_GRID_LEN][V_GRID_LEN] = {{2.98, 2.98, 2.98, 2.98, 2.98, 2.98}, 
                                                {2.98, 3.23, 3.23, 3.23, 3.23, 3.23}, 
                                                {3.23, 3.61, 3.61, 3.97, 3.70, 3.00}, 
                                                {3.23, 3.58, 3.87, 4.00,   4.,   4.}, 
                                                {3.23,   4.,   4.,  3.2, 2.41, 2.53}, 
                                                { 3.8,  3.8,   4.,   4.,   4.,   4.}, 
                                                { 2.4,  3.8,  3.8, 3.21, 2.91, 3.61}, 
                                                { .75,  .75,  .75, 1.46, 1.46, 1.46}, 
                                                {  0.,   0.,   0.,   0.,   0.,   0.}, 
                                                { 1.4,  1.4,  1.4,   0.,   0.,   0.}};

   double T1_GRID[YE_GRID_LEN][V_GRID_LEN] = {{.205, .205, .100, .080, .028, .028}, 
                                              {.321, .103, .066, .012, .012, .012}, 
                                              {.186, .122, .052, .026, .020, .020}, 
                                              {.186, .093, .047, .037, .021, .022}, 
                                              {.186, .131, .087, .052, .043, .015}, 
                                              { .14, .123, .089, .055, .045, .031}, 
                                              {.264,   .1,  .07, .055, .042, .033}, 
                                              {  0.,   0.,   0., .015, .013, .013}, 
                                              {  0.,   0.,   0.,   0.,   0.,   0.}, 
                                              {.051, .027, .011,   0.,   0.,   0.}}; 

   double SIG1_GRID[YE_GRID_LEN][V_GRID_LEN]={{ .16,  .18,  .10,  .20, .062, .062}, 
                                              {.225, .061, .129, .085, .085, .085}, 
                                              {.095, .095, .034, .020, .010, .010}, 
                                              {.095, .055,  .02,  .05,  .05,  .02}, 
                                              {.095,  .09, .051,  .02,  .02, .015}, 
                                              {.055, .067, .053, .032, .032, .024}, 
                                              {.075, .044,  .03,  .02,  .02, .014}, 
                                              {  .3,   .3,   .3, .004,.0015,.0015}, 
                                              {  .1,   .1,   .1,   .1,   .1,   .1}, 
                                              { .01,.0033, .002,   .1,   .1,   .1}};

   double C1_GRID[YE_GRID_LEN][V_GRID_LEN] = {{    0.,     0.,     0.,     0.,     0.,     0.}, 
                                              {    0.,     0.,     0.,     0.,     0.,     0.}, 
                                              {    0.,     0.,     0.,     0.,     0.,     0.}, 
                                              {28.901, 28.901, 28.138, 27.480, 27.480, 27.480}, 
                                              {27.909, 27.909, 27.086, 29.888, 30.575, 29.012}, 
                                              {27.539, 27.539, 26.525, 26.940, 29.289,     0.}, 
                                              {    0., 26.938, 28.036, 28.036, 27.447, 26.663}, 
                                              {35.576, 35.576, 35.576, 35.459, 35.459, 34.996}, 
                                              {20.723, 20.723, 20.906, 20.906, 20.906, 20.906}, 
                                              {25.786, 28.689, 29.012, 29.012, 29.710, 29.710}};

   double TAU1_GRID[YE_GRID_LEN][V_GRID_LEN] = {{   1.,   1.,   1.,   1.,   1.,   1.}, 
                                                {   1.,   1.,   1.,   1.,   1.,   1.}, 
                                                {   1.,   1.,   1.,   1.,   1.,   1.}, 
                                                { 4.15, 4.15, 5.36, 4.75, 4.75, .432}, 
                                                { 14.9, 11.2, 13.2, 5.96, 3.02, 5.35}, 
                                                { 12.2, 12.2, 17.2, 1.03, .613,.0086}, 
                                                {.0086, 13.8, 11.4, 14.3, 13.3, 13.3}, 
                                                { .038, .038, .038, .038, .029,  .06}, 
                                                { 1.29, 1.29, 1.29, 1.29, 1.29, 1.29}, 
                                                { 12.2, 1.38, 1.03, 1.03, 1.03, 1.03}};

   double C2_GRID[YE_GRID_LEN][V_GRID_LEN] = {{    0.,     0.,     0.,     0.,     0.,     0.}, 
                                              {    0.,     0.,     0.,     0.,     0.,     0.}, 
                                              {    0.,     0.,     0.,     0.,     0.,     0.}, 
                                              {23.273, 23.445, 23.157, 22.803, 22.803, 22.803}, 
                                              {21.416, 22.333, 22.481, 23.075, 24.407, 23.673}, 
                                              {    0.,     0.,     0., 19.802, 22.012, 21.045}, 
                                              {    0., 14.138, 18.792, 19.114, 23.810, 19.163}, 
                                              {28.782, 28.782, 28.324, 30.391, 30.073, 29.428}, 
                                              {16.760, 16.811, 16.811, 16.811, 16.811, 16.811}, 
                                              {23.901, 23.901, 23.901, 23.901, 23.901, 23.901}}; 
   
   double TAU2_GRID[YE_GRID_LEN][V_GRID_LEN] = {{  1.,   1.,   1.,   1.,   1.,   1.}, 
                                                {  1.,   1.,   1.,   1.,   1.,   1.}, 
                                                {  1.,   1.,   1.,   1.,   1.,   1.}, 
                                                {5.18, 4.14, 3.71, 3.37, 3.37, 3.37}, 
                                                {5.53, 6.91, 9.15, 5.01, 4.49, 4.49}, 
                                                {5.18, 5.18, 5.18, 34.7, 8.38, 22.6}, 
                                                {.864, 4.49, 95.0, 95.0, 0.95, 146.}, 
                                                {.043, .043, .086, .025, .034, .064}, 
                                                {3.45, 3.45, 3.45, 3.45, 3.45, 3.45}, 
                                                {25.9, 25.9, 25.9, 25.9, 25.9, 25.9}}; 

/*
  !************************************************************************
  !  Heating rates restricted to {v, ye}-grid.                            *
  !************************************************************************
*/  
 void heating_rate_grid(std::vector<body*>& bodies){
 void heating_rate_grid(int iv, int jye, double t){
      
   double e0, alp, t0, sig, alp1, t1, sig1, C1, C2, tau1, tau2;
   double a, b;
   double oneoverpi = 1./M_PI;
  
     e0=     E0_GRID[iv][jye];
     alp=   ALP_GRID[iv][jye];
     t0=     T0_GRID[iv][jye];
     sig=   SIG_GRID[iv][jye];
     alp1= ALP1_GRID[iv][jye];
     t1=     T1_GRID[iv][jye];
     sig1= SIG1_GRID[iv][jye];
     C1=     C1_GRID[iv][jye];
     tau1= TAU1_GRID[iv][jye];
     C2=     C2_GRID[iv][jye];
     tau2= TAU2_GRID[iv][jye];

     a= .5 - oneoverpi*atan((t - t0)/sig);
     b= .5d0 + oneoverpi*atan((t - t1)/sig1);
     double h= e0*1e18*(pow(a,alp) * pow(b,alp1)) 
               + exp(C1 - t/tau1*1e-3) + exp(C2 - t/tau2*1e-5);
 }


/*
  !************************************************************************
  !  Heating rates for arbitrary {v, Ye} at arbitrary time t[s]           *
  !************************************************************************
*/
 void heating_rate(double v, double ye, double t){

  /*
   v : ejecta expansion velocity [c]
   ye : initial electron fraction
   t : time [s]
  */

  int i1, i2, j1, j2;
  double v1, v2, y1, y2, fv, fy, f11, f12, f21, f22;
  double a, b;
  double oneoverpi = 1./M_PI;

  // Find index for v
     for (i1 = 0; i1 < V_GRID_LEN-1; ++i1){
       if (v < V_GRID[i1+1]){
         break;    
       } else if (i1 == 0 || i1 == V_GRID_LEN) {
          std::cout<<"ERROR : v outside the grid"<<std::endl;
          assert(false); 
       }
     }
     i2= i1 + 1;

  // Find index for ye
     for (j1 = 0; j1 < YE_GRID_LEN-1; ++j1){
       if (v < YE_GRID[j1+1]){
         break;    
       } else if (j1 == 0 || j1 == YE_GRID_LEN) {
          std::cout<<"ERROR : Ye outside the grid"<<std::endl;
          assert(false); 
       }
     }
     j2= j1 + 1;

     v1= V_GRID[i1];
     v2= V_GRID[i2];
     fv= (v - v1)/(v2 - v1);
     
     y1= YE_GRID[j1];
     y2= YE_GRID[j2];
     fy= (ye - y1)/(y2 - y1);

     f11= (1d0 - fv)*(1d0 - fy);
     f12= (1d0 - fv)*fy;
     f21= fv*(1d0 - fy);
     f22= fv*fy;

     e0=   f11*E0_GRID[i1][j1] + f12*E0_GRID[i1][j2] 
         + f21*E0_GRID[i2][j1] + f22*E0_GRID[i2][j2];

     alp=  f11*ALP_GRID[i1][j1] + f12*ALP_GRID[i1][j2] 
         + f21*ALP_GRID[i2][j1] + f22*ALP_GRID[i2][j2];

     t0=   f11*T0_GRID[i1][j1] + f12*T0_GRID[i1][j2] 
         + f21*T0_GRID[i2][j1] + f22*T0_GRID[i2][j2];

     sig=  f11*SIG_GRID[i1][j1] + f12*SIG_GRID[i1][j2] 
         + f21*SIG_GRID[i2][j1] + f22*SIG_GRID[i2][j2];

     alp1= f11*ALP1_GRID[i1][j1] + f12*ALP1_GRID[i1][j2] 
         + f21*ALP1_GRID[i2][j1] + f22*ALP1_GRID[i2][j2];

     t1=   f11*T1_GRID[i1][j1] + f12*T1_GRID[i1][j2] 
         + f21*T1_GRID[i2][j1] + f22*T1_GRID[i2][j2];

     sig1= f11*SIG1_GRID[i1][j1] + f12*SIG1_GRID[i1][j2] 
         + f21*SIG1_GRID[i2][j1] + f22*SIG1_GRID[i2][j2];

     C1=   f11*C1_GRID[i1][j1] + f12*C1_GRID[i1][j2] 
         + f21*C1_GRID[i2][j1] + f22*C1_GRID[i2][j2];

     tau1= f11*TAU1_GRID[i1][j1] + f12*TAU1_GRID[i1][j2] 
         + f21*TAU1_GRID[i2][j1] + f22*TAU1_GRID[i2][j2];

     C2=   f11*C2_GRID[i1][j1] + f12*C2_GRID[i1][j2] 
         + f21*C2_GRID[i2][j1] + f22*C2_GRID[i2][j2];

     tau2= f11*TAU2_GRID[i1][j1] + f12*TAU2_GRID[i1][j2] 
         + f21*TAU2_GRID[i2][j1] + f22*TAU2_GRID[i2][j2];

     a= .5 - oneoverpi*atan((t - t0)/sig);
     b= .5 + oneoverpi*atan((t - t1)/sig1);
     h= e0*1e18*(pow(a,alp) * pow(b,alp1)) 
      + exp(C1 - t/tau1*1e-3) + exp(C2 - t/tau2*1e-5);
 }



#if 0

/*
  !************************************************************************
  !  Output heating rates for arbitrary {v, Ye} in SuperNu format         *
  !************************************************************************
*/ 

  SUBROUTINE output_supernu(v, ye)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN):: v  ! ejecta expansion velocity [c]
  DOUBLE PRECISION, INTENT(IN):: ye ! initial electron fraction
  !
  INTEGER, PARAMETER :: nmax = 1000
  DOUBLE PRECISION, PARAMETER:: &
    a=  0.05d0,      & !  5% alphas
    b=  0.20d0,      & ! 20% betas
    nu= 0.35d0,      & ! 35% neutrinos
    rad= a + b + nu, & ! gammas
    gam= 1d0 - rad,  & ! gammas
    tmin= 1d-3,      & ! [s] initial time
    tmax= 86400*100    ! [s] final time (100 days)
  INTEGER :: n
  DOUBLE PRECISION:: t, dtfac, hr

  PRINT ('("# Heating rates from Stephan Rosswog (WinNet),",I2,"% alpha'// &
         's, ",I2,"% betas, ",I2,"% nu")'),INT(a*100),INT(b*100),INT(nu*100)
  PRINT ('("# 1st block is thermalization parameters, 2nd block is heat'// &
         'ing rates")')
  PRINT *
  PRINT ('("# 1st line is thermalization option for columns 4-8 in the '// &
         'heating rates")')
  PRINT ('("# 2nd line is thermalization parameter (use depends on opti'// &
         'on in 1st line)")')
  PRINT ('("           1             1             0             1     '// &
         '        1")')
  PRINT ('("1.200000E-11  1.300000E-11  1.000000E+00  1.300000E-11  2.0'// &
         '00000E-12")')
  PRINT *
  PRINT ('("# ",I4)'), nmax
  PRINT ('("# 1:time[s], 2:total nuclear heating rate[erg/(g*s)] 3:tota'// &
         'l radiation[erg/g*s]")')
  PRINT ('("# 4:alpha 5:beta 6:gamma 7:electrons 8:fission products")')

  t= tmin
  dtfac= exp((log(tmax) - log(tmin))/(dble(nmax) - 1))
  DO n=1,1000
     hr= heating_rate(v, ye, t)
     PRINT '(12(ES12.5,1X))', t, hr, hr*rad, hr*a, hr*b, hr*gam, 0d0, 0d0
     t= t*dtfac
  ENDDO

  END SUBROUTINE output_supernu

ENDMODULE hratelib

! simple test and usage example
PROGRAM heating_rate_test
USE hratelib
IMPLICIT NONE
INTEGER :: n
INTEGER, PARAMETER :: nmax = 1000
DOUBLE PRECISION:: v, ye, t, dtfac
DOUBLE PRECISION, PARAMETER:: tmin= 1d-3       ! [s]
DOUBLE PRECISION, PARAMETER:: tmax= 86400*100  ! 100 days [s]
CHARACTER(len=32):: arg

   ! help message
   IF (iargc().NE.4) STOP "Usage example: ./hrate.x -v 0.2 -ye 0.25"

   argparse: DO n = 1, iargc()
      CALL getarg(n, arg)
      IF (TRIM(arg).EQ."-v") THEN
         CALL getarg(n+1, arg)
         READ(arg, *) v
      ENDIF
      IF (TRIM(arg).EQ."-ye") THEN
         CALL getarg(n+1, arg)
         READ(arg, *) ye
      ENDIF
   ENDDO argparse

   CALL output_supernu(v, ye)
   !! ! output time vs heating rate only
   !!
   !! PRINT '(A)', "# 1:time[s] 2:heating rate[erg/g/s]"
   !! t= tmin
   !! dtfac= exp((log(tmax) - log(tmin))/(dble(nmax) - 1))
   !! DO n=1,1000
   !!    PRINT '(2(ES12.5,1X))', t, heating_rate(v, ye, t)
   !!    t= t*dtfac
   !! ENDDO
#endif

#endif
