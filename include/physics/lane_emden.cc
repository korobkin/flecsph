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

#include <iostream>
#include <cmath>
#include <vector>




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
double A_wd  = pi*m_e*pow(clite,2)/(3.0*pow(lam_e,3));
double B_wd  = 8.0*pi/(3.0*N_A*pow(lam_e,3));

// --- initialization
//rho_c = (1.2*Msun)/(4*pi/3*(5000*km)**3)  // [g/cm3] central density
double rho_c = 1.822425e6; // [g/cm3] central density

// int Nr = 2000;          // number of intervals in radial direction
double Y_e = 0.5;          // electron fraction



// --- equation of state
double eos_pressure(double rho){
	//return K*pow(rho,gam) # polytropic
	double x = pow(rho*Y_e/B_wd, 1./3.);
	return A_wd*(x*sqrt(1 + pow(x,2))*(2*pow(x,2) - 3) + 3*asinh(x));
}
  

double eos_dPdrho(double rho){
	//return gam*eos_pressure(rho)/rho
	double x = pow(rho*Y_e/B_wd, 1./3.);
	return 8.0*A_wd*Y_e*pow(x,2)/(3*B_wd*sqrt(1 + pow(x,2)));
}
  

double eos_eint(double rho){
	//return K*pow(rho,(gam-2))/(gam-1); // polytropic
	double x = pow(rho*Y_e/B_wd, 1./3.);
	double x2 = x*x;
	double x3 = x*x2;
	return A_wd/rho*(8*x3*(sqrt(1 + x2) - 1)
                 - (x*(2*x2 - 3.)*sqrt(1 + x2) + 3.*asinh(x)));
}

// rho = rho_c * theta**n
double n = 1.0/(rho_c / eos_pressure(rho_c) * eos_dPdrho(rho_c) - 1.0);

// pseudo polytropic EOS for first step 
double K = eos_pressure(rho_c) / pow(rho_c,(1.0 + 1.0/n));
double gam = 1 + 1.0/n;

// useful constant for the first step
double alpha = 4*pi*ggrav / (K*(n+1)*pow(rho_c,(1.0/n)));


// ds/dtheta (s = r**2)
double dsdtheta(double m, double s, double theta){
    double rho = rho_c * pow(theta, n);
    return -2*n*pow(s,1.5)/(ggrav * m * theta) * eos_dPdrho(rho);
}

// dm/dtheta
double dmdtheta(double m, double s, double theta){
    double rho = rho_c * pow(theta, n);
    return -4*pi*n*rho_c*pow(theta, n-1)*pow(s,2) / (ggrav * m) * eos_dPdrho(rho);
}

// one step for RK4 integration for the Lane-Emden solver
std::pair<double, double> lane_emden_RK4(double m, double s, double theta, double dtheta){
	double km_1 =  dmdtheta(m,s,theta);
    double ks_1 =  dsdtheta(m,s,theta);
    
    double km_2 =  dmdtheta(m + dtheta/2*km_1, s + dtheta/2*ks_1, theta + dtheta/2);
    double ks_2 =  dsdtheta(m + dtheta/2*km_1, s + dtheta/2*ks_1, theta + dtheta/2);
    
    double km_3 =  dmdtheta(m + dtheta/2*km_2, s + dtheta/2*ks_2, theta + dtheta/2);
    double ks_3 =  dsdtheta(m + dtheta/2*km_2, s + dtheta/2*ks_2, theta + dtheta/2);
    
    double km_4 =  dmdtheta(m + dtheta*km_3, s + dtheta*ks_3, theta + dtheta);
    double ks_4 =  dsdtheta(m + dtheta*km_3, s + dtheta*ks_3, theta + dtheta);
    
    double m_ret = m + dtheta/6*(km_1 + 2*km_2 + 2*km_3 + km_4);
    double s_ret = s + dtheta/6*(ks_1 + 2*ks_2 + 2*ks_3 + ks_4);
	return std::make_pair(m_ret, s_ret);
}

// Main function for solving Lane-Emden equation and printing out the result for storing
void lane_emden_solver(int Nr){
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


	
	// RK4 for integration of the two ODEs
	for(int i = 0; i < Nr-1; i++){
		theta_cur = theta_arr[i];
	    std::pair<double, double> lane_emden_result = lane_emden_RK4(m_arr[i], s_arr[i], theta_arr[i], theta_step);
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
	printf ("#  - central pressure:  %12.12e [dynes/cm^2]\n", eos_pressure(rho_c));
	printf ("#\n");
	printf ("# Equation of state: zero-temperature WD\n");
	// printf ("#  - K = %12.5e" % K)
	// printf ("#  - \\Gamma = %12.5e" % gam)

	// normalization
	double rho_norm = M_star/pow(R_star,3);
	double drhodr_norm = M_star/pow(R_star,4);

	//printing the density profile for the .dat file
	printf("%19.12e %19.12e %19.12e %19.12e\n", 0.0, rho_c/rho_norm, 0.0, 0.0);
	for(int i = 0; i < Nr - 1; i++){
		double m = m_arr[i];
		double r = sqrt(s_arr[i]);
		double rho = rho_c * pow(theta_arr[i],n);
		double drhodr = -ggrav*m*rho/(r*r * eos_dPdrho(rho));
		
		printf("%19.12e %19.12e %19.12e %19.12e\n", r/R_star, rho/rho_norm, m/M_star, drhodr/drhodr_norm);
	}
	printf("%19.12e %19.12e %19.12e %19.12e\n", 1.0, 0.0, 1.0, 0.0);
}


int main(){
	int Nr = 10000;
	lane_emden_solver(Nr);
	
	return 0;
}
