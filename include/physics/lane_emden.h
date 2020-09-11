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

namespace lane_emden{
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
			particle.setDensity(rho);
			eos::compute_pressure(particle);
			return particle.getPressure();
		}


		double eos_dPdrho(body& particle, double rho){
			particle.setDensity(rho);
			eos::compute_soundspeed(particle);
			double cs = particle.getSoundspeed();
			return cs*cs*rho;
		}


		double eos_eint(body& particle, double rho){
			particle.setDensity(rho);
			eos::compute_internal_energy(particle);
			return particle.getInternalenergy();
		}
		
		// ds/dtheta (s = r**2)
		double dsdtheta(double m, double s, double theta, double rho_c, body& pt){
			double n = 1.0/(rho_c / eos_pressure(pt, rho_c) * eos_dPdrho(pt, rho_c) - 1.0);
		    double rho = rho_c * pow(theta, n);
		    return -2*n*pow(s,1.5)/(ggrav * m * theta) * eos_dPdrho(pt, rho);
		}

		// dm/dtheta
		double dmdtheta(double m, double s, double theta, double rho_c, body& pt){
			double n = 1.0/(rho_c / eos_pressure(pt, rho_c) * eos_dPdrho(pt, rho_c) - 1.0);
		    double rho = rho_c * pow(theta, n);
		    return -4*pi*n*rho_c*pow(theta, n-1)*pow(s,2) / (ggrav * m) * eos_dPdrho(pt, rho);
		}

		// one step for RK4 integration for the Lane-Emden solver
		std::pair<double, double> lane_emden_RK4(double m, double s, double theta, double dtheta, double rho_c, body& pt){
			double km_1 =  dmdtheta(m,s,theta,rho_c, pt);
		    double ks_1 =  dsdtheta(m,s,theta,rho_c, pt);
		    
		    double km_2 =  dmdtheta(m + dtheta/2*km_1, s + dtheta/2*ks_1, theta + dtheta/2,rho_c, pt);
		    double ks_2 =  dsdtheta(m + dtheta/2*km_1, s + dtheta/2*ks_1, theta + dtheta/2,rho_c, pt);
		    
		    double km_3 =  dmdtheta(m + dtheta/2*km_2, s + dtheta/2*ks_2, theta + dtheta/2,rho_c, pt);
		    double ks_3 =  dsdtheta(m + dtheta/2*km_2, s + dtheta/2*ks_2, theta + dtheta/2,rho_c, pt);
		    
		    double km_4 =  dmdtheta(m + dtheta*km_3, s + dtheta*ks_3, theta + dtheta,rho_c, pt);
		    double ks_4 =  dsdtheta(m + dtheta*km_3, s + dtheta*ks_3, theta + dtheta,rho_c, pt);
		    
		    double m_ret = m + dtheta/6*(km_1 + 2*km_2 + 2*km_3 + km_4);
		    double s_ret = s + dtheta/6*(ks_1 + 2*ks_2 + 2*ks_3 + ks_4);
			return std::make_pair(m_ret, s_ret);
		}

		

		std::vector< std::vector<double> > lane_emden(double rho_c,
													  int Nr){
													  // double (*eos_pressure)(double),
													  // double (*eos_dPdrho)(double),
													  // double (*eos_eint)(double) ){
													  // std::function<double(double)> eos_pressure,
													  // std::function<double(double)> eos_dPdrho,
													  // std::function<double(double)> eos_eint){
			eos::select();

			body pt0;
			eos::eos_init(pt0);

			// rho = rho_c * theta**n
			double n = 1.0/(rho_c / eos_pressure(pt0, rho_c) * eos_dPdrho(pt0, rho_c) - 1.0);

			// pseudo polytropic EOS for first step 
			double K = eos_pressure(pt0, rho_c) / pow(rho_c,(1.0 + 1.0/n));
			double gam = 1 + 1.0/n;

			// useful constant for the first step
			double alpha = 4*pi*ggrav / (K*(n+1)*pow(rho_c,(1.0/n)));

			// define return variables
			std::vector< std::vector<double> > ret2dvec;
			std::vector<double> rad_arr(Nr+1, 0);
			std::vector<double> rho_arr(Nr+1, 0);
			std::vector<double> mass_arr(Nr+1, 0);
			std::vector<double> drhodr_arr(Nr+1, 0);

			//start the solver
			double theta_min = 1e-7 / Nr;
			double theta_step = - (1.0 - theta_min) / (Nr-1);
			std::vector<double> theta_arr(Nr);
			for(int i = 0; i < Nr; i++) {
				theta_arr[i] = 1.0 + i * theta_step;
			}

			// first step is approximated with polytropic EOS with const rho = rho_c, which gives
			// dm = 4*pi/3*rho_c*dr**3, ds = -6.0/(alpha*rho_c) * theta_step;
			double s_init = -6.0/(alpha*rho_c) * theta_step;
			double m_init = 4.0*pi/3*pow(s_init,1.5) * rho_c;
			double theta_cur = 1.0;

			std::vector<double> s_arr(Nr);
			std::vector<double> m_arr(Nr);
			s_arr[0] = s_init;
			m_arr[0] = m_init;

			std::cout<<"rhoc: "<<rho_c<<"p_c"<<eos_pressure(pt0,rho_c)<<std::endl;
			
			// RK4 for integration of the two ODEs
			for(int i = 0; i < Nr-1; i++){
				theta_cur = theta_arr[i];
			    std::pair<double, double> lane_emden_result = lane_emden_RK4(m_arr[i], s_arr[i], theta_arr[i], theta_step, rho_c, pt0);
			    m_arr[i+1] = lane_emden_result.first;
			    s_arr[i+1] = lane_emden_result.second;
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

			rad_arr[Nr] = 1.0;
			rho_arr[0] = rho_c / rho_norm;
			mass_arr[Nr] = 1.0;


			for(int i = 0; i < Nr - 1; i++){
				double m = m_arr[i];
				double r = sqrt(s_arr[i]);
				double rho = rho_c * pow(theta_arr[i],n);
				double drhodr = -ggrav*m*rho/(r*r * eos_dPdrho(pt0,rho));
				mass_arr[i+1] = m / M_star;
				rad_arr[i+1] = r / R_star;
				rho_arr[i+1] = rho / rho_norm;
				drhodr_arr[i+1] = drhodr/drhodr_norm;
				// printf("%19.12e %19.12e %19.12e %19.12e\n", r/R_star, rho/rho_norm, m/M_star, drhodr/drhodr_norm);
			}



			ret2dvec.push_back(rad_arr);
			ret2dvec.push_back(rho_arr);
			ret2dvec.push_back(mass_arr);
			ret2dvec.push_back(drhodr_arr);
			return ret2dvec;
		}
}



#endif // lane_emden.h