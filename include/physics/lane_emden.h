/** This file contains a Lane-Emden solver for white dwarf (WD) or polytropic EOS (currently for WD)
 *  instead of evolving rho, m with drho/dr, dm/dr -> evolve s,m with ds/dtheta, dm/dtheta
 *  (s = r**2, rho = rho_c * theta**n)
 *  The inverse integration allows both ends of the integrated variable to be defined: theta \in [1.0, 0.0]
 *  The use of s = r**2 allows an analytical form for approximating the first step where r = 0.
 *
 *  The code shows convergence rate = 1.5.
 *  Use: $ g++ -std=c++11 lane_emden.cc -o <excutable filename>
 *       $ ./<excutable filename> > <.dat filename>
 *  (Johnny 09.02.2020)
 */
#ifndef LANE_EMDEN_H
#define LANE_EMDEN_H

#include <iostream>
#include <cmath>
#include <vector>

#include "eos.h"
#include "body.h"

namespace lane_emden {

// --- physics constants
double clite = 2.99792458e10 ;  // [cm/s]
double ggrav = 6.67408e-8    ;  // [cm3/(g s^2)]
double Msun  = 1.9891e+33    ;  // [g]
double km    = 1.0e5         ;  // [cm]
double hpl   = 6.62607015e-27;  // [erg*s]
double m_e   = 9.10938370e-28;  // [g]
double N_A   = 6.02214076e+23;  // [mol^-1]
double pi    = 3.141592653589;

double lam_e = hpl/m_e/clite;  // [cm]

// --- equation of state
double eos_pressure(body& particle, double rho){
  // std::cout<<"setting pressure, density is "<< rho <<std::endl;
  particle.setDensity(rho);
  eos::compute_pressure(particle);
  return particle.getPressure();
}

double eos_dPdrho(body& particle, double rho){
  particle.setDensity(rho);
  eos::compute_soundspeed(particle);
  double cs = particle.getSoundspeed();
  return cs*cs;
}

double eos_eint(body& particle, double rho){
  particle.setDensity(rho);
  eos::compute_internal_energy(particle);
  return particle.getInternalenergy();
}
		
// ds/dth (s = r**2, and th is defined from rho = rho_c\th^n)
std::pair<double, double>
dms_dth(double m, double s, double th, double rho_c, double n, body& pt){
  double rho = rho_c * pow(th, n);
  pt.setDensity(rho);
  eos::compute_soundspeed(pt);
  double cs = pt.getSoundspeed();
  double dPdrho_S = cs*cs;
  double dsdth = -2*n*sqrt(s*s*s)/(ggrav * m * th) * dPdrho_S;
  double dmdth = dsdth * 2*pi*sqrt(s)*rho;
  return {dmdth, dsdth};
}

// one step for RK4 integration for the Lane-Emden solver
std::pair<double, double> 
lane_emden_RK4(double m, double s,
    double th, double dth, double rho_c, double n, body& pt){
  double th12 = th + dth/2;
  double th2  = th + dth;
  auto [km_1, ks_1] = dms_dth(m,s,th,rho_c,n, pt);
  auto [km_2, ks_2] = dms_dth(m + dth/2*km_1,s + dth/2*ks_1,th12,rho_c,n,pt);
  auto [km_3, ks_3] = dms_dth(m + dth/2*km_2,s + dth/2*ks_2,th12,rho_c,n,pt);
  auto [km_4, ks_4] = dms_dth(m + dth*km_3,s + dth*ks_3,th2,rho_c,n,pt);
  
  double m_ret = m + dth/6*(km_1 + 2*km_2 + 2*km_3 + km_4);
  double s_ret = s + dth/6*(ks_1 + 2*ks_2 + 2*ks_3 + ks_4);

  return {m_ret, s_ret};
}

void solve(double rho_c, double p_c, int Nr,
    std::vector<double> & rad_arr,
    std::vector<double> & rho_arr,
    std::vector<double> & mass_arr,
    std::vector<double> & drhodr_arr) {

  // eos::select();

  body pt0;
  pt0.setDensity(rho_c);
  pt0.setPressure(p_c);
  eos::compute_entropy(pt0);
  eos::compute_soundspeed(pt0);
  double cs = pt0.getSoundspeed();
  double dPdrho_c = cs*cs;

  // rho = rho_c * theta**n
  double gam = rho_c/p_c*dPdrho_c;
  double n = 1./(gam - 1.);

  // pseudo polytropic EOS for first step 
  double K_c = p_c / pow(rho_c, gam);

  // useful constant for the first step
  double alpha = 4*pi*ggrav / (K_c*(n+1)*pow(rho_c,(1.0/n)));

  // allocate arrays
  rad_arr.resize(Nr);
  rho_arr.resize(Nr);
  mass_arr.resize(Nr);
  drhodr_arr.resize(Nr);

  //start the solver
  double theta_min = 1e-7 / (double)Nr;  // TODO: make 1e-7 into a parameter
  double theta_step = - (1.0 - theta_min)/(double)(Nr - 1);
  std::vector<double> theta_arr(Nr);
  for(int i = 0; i < Nr; i++) {
      theta_arr[i] = 1.0 + i*theta_step;
  }

  // first step is approximated with polytropic EOS with const rho = rho_c, which gives
  // dm = 4*pi/3*rho_c*dr**3, ds = -6.0/(alpha*rho_c) * theta_step;
  double s_init =-6./(alpha*rho_c) * theta_step;
  double m_init = 4.*pi/3*sqrt(CU(s_init)) * rho_c;
  double theta_cur = 1.;

  std::vector<double> s_arr(Nr);
  std::vector<double> m_arr(Nr);
  s_arr[0] = 0.;
  m_arr[0] = 0.;
  s_arr[1] = s_init;
  m_arr[1] = m_init;

  // RK4 for integration of the two ODEs
  for(int i = 1; i < Nr - 1; i++) {
    theta_cur = theta_arr[i];
    auto [first, second] = lane_emden_RK4(m_arr[i],s_arr[i], 
        theta_arr[i], theta_step, rho_c, n, pt0);
    m_arr[i+1] = first;
    s_arr[i+1] = second;
  }
  
  // Finally!
  double M_star = m_arr[Nr-1];
  double R_star = sqrt(s_arr[Nr-1]);

  //printing the header for the .dat file
  printf ("# Stellar parameters:\n");
  printf ("#  - mass:    %12.12e [g]\n", M_star);
  printf ("#  - radius:  %12.12e [cm]\n", R_star);
  printf ("#  - central density:  %12.12e [g/cm^3]\n", rho_c);
  printf ("#  - central pressure:  %12.12e [dynes/cm^2]\n", eos_pressure(pt0,rho_c));
  printf ("#\n");
  printf ("# Equation of state: zero-temperature WD\n");

  // normalization
  double rho_norm = M_star/pow(R_star,3);
  double drhodr_norm = M_star/pow(R_star,4);

  for(int i = 0; i < Nr; i++){
    double m = m_arr[i];
    double r = sqrt(s_arr[i]);
    double rho = rho_c * pow(theta_arr[i],n);
    double drhodr = -ggrav*m*rho/(r*r * eos_dPdrho(pt0,rho));
    mass_arr[i] = m / M_star;
    rad_arr[i] = r / R_star;
    rho_arr[i] = rho / rho_norm;
    drhodr_arr[i] = (i ? drhodr/drhodr_norm : 0.);
  }

  rad_arr[0] = 0.;
  rho_arr[0] = rho_c / rho_norm;
  mass_arr[0] = 0.;
  drhodr_arr[0] = 0.;

  rad_arr[Nr - 1] = 1.;
  rho_arr[Nr - 1] = 0.;
  mass_arr[Nr - 1] = 1.;
  drhodr_arr[Nr - 1] = 0.;

  // TODO: 1. add a parameter: string 'lane_emden_output_profile'
  //       2. only output from MPI rank 0
  //       2. if string is empty (zero length), do not output profile;
  //       3. if string is non-empty, assume it contains profile file name;
  //       4. attempt to create file with that name;
  //       5. if the file already exists, issue a warning and overwrite it;
  //       6. check that the file has been successfully created;
  //       7. output the profile data using format below
  for(int i = 0; i < Nr; i++){
    printf("%19.12e %19.12e %19.12e %19.12e\n", 
        rad_arr[i], rho_arr[i], mass_arr[i], drhodr_arr[i]);
  }
} // solve(..)

} // namespace lane_emden



#endif // lane_emden.h
