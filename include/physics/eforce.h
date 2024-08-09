/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

/*~--------------------------------------------------------------------------~*
 *
 * /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file eforce.h
 * @brief Namespace for the choice of external force and external potential
 */

#ifndef _eforce_h_
#define _eforce_h_

#include <variant>
#include <boost/algorithm/string.hpp>

#include "density_profiles.h"
#include "params.h"
#include "tree.h"
#include "kernels.h"
#include "math.h"

#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))

namespace external_force {


/*-----------------------------------------------------------------------------*
 * class force_base
 * @brief implementation interface for external forces
 *-----------------------------------------------------------------------------*/
template<class Derived>
struct force_base {
  /**
  * @brief      interface of potential
  *
  * @param[in]  p   coordinates to evaluate potential
  */
  inline double potential(const point_t& p) const
  {
    return static_cast<Derived>(*this).potential(p);
  }

  /**
  * @brief      interface of acceleration
  *
  * @param[in]  b   vector to evaluate acceleration
  */
  inline point_t acceleration(const body& b) const
  {
    return static_cast<Derived>(*this).acceleration(b);
  }
  private:
    force_base() = default;
    friend Derived;
};

/**
 * @brief      1D walls: steep power-law-like potentials
 */
template<auto I = 0>
struct force_square_well : public force_base<force_square_well<I>> {
  const double box[3] = {.5*param::box_length,
                         .5*param::box_width,
                         .5*param::box_height},
               pw_n = param::extforce_wall_powerindex,
               pw_a = param::extforce_wall_steepness;

  inline double potential(const point_t& rp) const {
    return pw_a*(((rp[I] <-box[I]) ? pow(-rp[I] - box[I], pw_n) : 0.0) +
                 ((rp[I] > box[I]) ? pow( rp[I] - box[I], pw_n) : 0.0));
  }


  inline point_t acceleration(const body & b ) const {
    point_t a = 0.0;
    point_t rp = b.coordinates();
    a[I] = pw_a*pw_n*(((rp[I] <-box[I])? pow(-rp[I] - box[I], pw_n - 1):0) -
                      ((rp[I] > box[I])? pow( rp[I] - box[I], pw_n - 1):0));
    return a;
   }
};


/**
 * @brief      Round or spherical boundary wall
 */
struct force_spherical_wall : public force_base<force_spherical_wall> {

  inline double potential(const point_t & rp) const {
    const double pw_n = param::extforce_wall_powerindex,
                 pw_a = param::extforce_wall_steepness,
                 R_sp = param::sphere_radius;
    double r = flecsi::magnitude(rp);
    return (r > R_sp) ? (pw_a*pow(r - R_sp, pw_n)) : 0.0;
  }


  inline point_t acceleration(const body & particle) const{
    point_t a = 0.0;
    point_t rp = particle.coordinates();
    const double pw_n = param::extforce_wall_powerindex,
                 pw_a = param::extforce_wall_steepness,
                 R_sp = param::sphere_radius;
    double r = flecsi::magnitude(rp);
    if(r > R_sp) {
      const double ar = pw_a*pw_n*pow(r - R_sp, pw_n - 1);
      for(unsigned short i = 0; i < gdimension; ++i)
        a[i] = -rp[i] / r * ar;
    }
    return a;
  }
};


/**
 * @brief    External force support for a spherically-symmetric
 *           density profile (from density_profiles.h)
 */
struct force_spherical_density_support :
public force_base<force_spherical_density_support> {

  force_spherical_wall _fsw;

  inline double potential(const point_t & rp) const {
    using namespace param;
    const double K0 = pressure_initial 
                    / pow(rho_initial, poly_gamma),
               rho0 = density_profiles::spherical_density_profile(0.),
                  x = flecsi::magnitude(rp) / sphere_radius;
    double rho = rho_initial / rho0
               * density_profiles::spherical_density_profile(x);
    double phi = (rho > 0)
        ? (-K0*poly_gamma*pow(rho, poly_gamma - 1.)/(poly_gamma - 1.))
        : 0;
    return phi + _fsw.potential(rp);
  }


  inline point_t acceleration(const body & particle) const {
    using namespace param;
    point_t a = 0.0;
    point_t rp = particle.coordinates();
    const double K0 = pressure_initial 
                    / pow(rho_initial, poly_gamma),
               rho0 = density_profiles::spherical_density_profile(0.),
                  r = flecsi::magnitude(rp),
                  x = r / sphere_radius;
    if(x > 1e-12) {
      double rho = rho_initial / rho0
                 * density_profiles::spherical_density_profile(x);
      double drhodr = rho_initial / (rho0 * sphere_radius)
                    * density_profiles::spherical_drho_dr(x);
      double a_r = (rho > 0)
          ? (K0*poly_gamma*pow(rho, poly_gamma - 2) * drhodr)
          : 0;
      for(short int i = 0; i < gdimension; ++i)
        a[i] = a_r * rp[i] / r;
    }
    return a + _fsw.acceleration(particle);
  }

};


/**
 * @brief    External force support for an arbitrary 3D density
 *           (from density_profiles.h)
 */
struct force_ndim_density_support :
public force_base<force_ndim_density_support> {

  force_spherical_wall _fsw;

  inline double potential(const point_t & rp) const {
    using namespace param;
    const double K0 = pressure_initial 
                    / pow(rho_initial, poly_gamma),
                rho = density_profiles::density_ndim(rp);
    double phi = (rho > 0)
        ? (-K0*poly_gamma*pow(rho, poly_gamma - 1.)/(poly_gamma - 1.))
        : 0;
    return phi + _fsw.potential(rp);
  }

  inline point_t acceleration(const body & particle) const {
    using namespace param;
    point_t rp = particle.coordinates();
    const double  K0 = pressure_initial 
                     / pow(rho_initial, poly_gamma),
                 rho = density_profiles::density_ndim(rp);
    point_t grad_rho = density_profiles::grad_density_ndim(rp);
    point_t a = (rho > 0) 
              ? (K0*poly_gamma*pow(rho, poly_gamma - 2) * grad_rho)
              : 0;
    return a + _fsw.acceleration(particle);
  }
};

/**
 * @brief    Add uniform constant gravity acceleration
 * 	         in y-direction (or x-direction if number of
 * 	         dimensions == 1)
 */
struct force_gravity : public force_base<force_gravity> {

  inline double potential(const point_t & rp) const {
    double height = rp[0];
    if(gdimension > 1)
      height = rp[1];
    return height * param::gravity_acceleration_constant;
  }


  inline point_t acceleration(const body & particle) const {
    point_t acc = 0.0;
    if(gdimension > 1)
      acc[1] = -param::gravity_acceleration_constant; // negative y-direction
    else
      acc[0] = -param::gravity_acceleration_constant; // negative x-direction
    return acc;
  }

};


/**
 * @brief    Gravitational field of a point mass at the origin, with a softening
 */
struct force_central_mass : public force_base<force_central_mass> {

  inline double potential(const point_t & rp) const {
    const double G = param::gravitational_constant;
    const double M = param::extforce_central_mass;
    const double eps = param::extforce_mass_softening_radius;
    double r = rp[0]*rp[0];
    for(unsigned short i = 1; i < gdimension; ++i)
      r += rp[i]*rp[i];
    r = sqrt(r + eps*eps);
    return -G*M/r;
  }


  inline point_t acceleration(const body & particle) const {
    const double G = param::gravitational_constant;
    const double M = param::extforce_central_mass;
    const double eps = param::extforce_mass_softening_radius;
    point_t rp = particle.coordinates(); 
    double r = rp[0]*rp[0];
    for(unsigned short i = 1; i < gdimension; ++i)
      r += rp[i]*rp[i];
    r = sqrt(r + eps*eps);
    return -G*M/(r*r*r)*rp;
  }

};


/**
 * @brief      2D airfoil in a wind tunnel
 *
 * The airfoil profile is centered at the anchor, tilted
 * at an angle to the flow. The shape of the airfoil can be described by the
 * following three parameters:
 *  - airfoil_size:           airfoil horizontal extent;
 *  - airfoil_thickness:      how thick is it;
 *  - airfoil_camber:         maximum deviation of camber line from the chord.
 *
 * Airfoil is positioned and rotated relative to its rear tip:
 *  - airfoil_anchor_x:       the x-coordinate of the anchor;
 *  - airfoil_anchor_y:       the y-coordinate of the anchor;
 *  - airfoil_attack_angle:   angle of attack - rotation from initial position
 *                            which is parallel to the x-axis.
 *
 * @param      particle  The particle being accelerated
 */
struct force_airfoil : public force_base<force_airfoil> {
  const double alpha = param::airfoil_attack_angle * M_PI / 180.0;

  inline double potential(const point_t & rp) const {
    using namespace param;
    double phi = 0.0;
    assert(gdimension > 1);

    static const double alpha = airfoil_attack_angle * M_PI / 180.0,
                        pw_n = extforce_wall_powerindex,
                        pw_a = extforce_wall_steepness;
    const double x1 = rp[0] - airfoil_anchor_x, y1 = rp[1] - airfoil_anchor_y;
    const double x = x1 * cos(alpha) + y1 * sin(alpha),
                 y = -x1 * sin(alpha) + y1 * cos(alpha);

    bool inside_bounding_box = std::abs(y) < 5.0 * airfoil_thickness &&
                               x > -airfoil_size * 0.02 &&
                               x < airfoil_size * 1.02;
    double upper_surface =
      airfoil_thickness * x * sqrt(airfoil_size * airfoil_size - x * x);
    double camber_line = airfoil_camber * sin(M_PI * x / 2.);
    double aux = SQ(upper_surface) - SQ(y - camber_line) + 0.002;
    if(inside_bounding_box && aux > 0.0)
      phi = pw_a * pow(aux, pw_n);
    return phi;
  }


  inline point_t acceleration(const body & particle) const {
    using namespace param;
    point_t a = 0.0;
    assert(gdimension > 1);

    point_t rp = particle.coordinates();
    const double x1 = rp[0] - airfoil_anchor_x,
                 y1 = rp[1] - airfoil_anchor_y,
                 alpha = airfoil_attack_angle * M_PI / 180.0,
                 pw_n = extforce_wall_powerindex,
                 pw_a = extforce_wall_steepness;
    const double x =  x1*cos(alpha) + y1*sin(alpha),
                 y = -x1*sin(alpha) + y1*cos(alpha);

    bool inside_bounding_box = std::abs(y) < 5.0 * airfoil_thickness &&
                               x > -airfoil_size * 0.02 &&
                               x <  airfoil_size * 1.02;
    double upper_surface =
      airfoil_thickness * x * sqrt(airfoil_size * airfoil_size - x * x);
    double camber_line = airfoil_camber * sin(M_PI * x / 2.);
    double phi = SQ(upper_surface) - SQ(y - camber_line) + 0.002;
    if(inside_bounding_box && phi > 0.0) {
      double a0, a1;
      a0 = pw_n * pw_a * pow(phi, pw_n - 1) *
           (2. * (y - camber_line) *
               (-airfoil_camber * M_PI / 2. * cos(M_PI / 2. * x)) -
             airfoil_thickness * airfoil_thickness * 2 * x *
               (airfoil_size * airfoil_size - 2 * x * x));
      a1 = pw_n * pw_a * pow(phi, pw_n - 1) * 2. * (y - camber_line);
      a[0] = a0 * cos(alpha) - a1 * sin(alpha);
      a[1] = a0 * sin(alpha) + a1 * cos(alpha);
    }
  }
};

struct
force_orbit : public force_base<force_orbit> {

  const double grav = param::gravitational_constant,
               a_sp = param::orbital_separation,
               m_1 = param::mass_primary_star,
               m_2 = param::mass_secondary_star,
               m_t = m_2 + m_1;

  inline double potential(const point_t & rp) const {
    using namespace param;
    assert(gdimension == 3); // TODO: generalize to 2D!
    double term1 = -grav*m_2/sqrt(SQ(rp[0] - a_sp) + SQ(rp[1]) + SQ(rp[2]));
    double term2 = -0.5*grav*m_t/CU(a_sp)
                 * (SQ(rp[0] - a_sp*m_2/m_t) + SQ(rp[1]));
    return term1 + term2;
  }


  inline point_t acceleration(const body & particle) const {
    using namespace param;
    assert(gdimension == 3); // TODO: generalize to 2D!
    point_t rp = particle.coordinates();
    point_t acc = 0.0;

    double temp = SQ(rp[0] - a_sp) + SQ(rp[1]) + SQ(rp[2]);
    temp = sqrt(CU(temp));
    double term1 = -grav * m_2 / temp;
    acc[0] += term1 * (rp[0] - a_sp);
    acc[1] += term1 * rp[1];
    acc[2] += term1 * rp[2];

    double term2 = grav * m_t / CU(a_sp);
    acc[0] += term2 * (rp[0] - a_sp * m_2 / m_t); // x-direction
    acc[1] += term2 * rp[1];
    return acc;
  }

};

/**
 * @brief      Constant potential shift: used for debugging
 */
struct force_poison : public force_base<force_poison> {

  inline double
  potential(const point_t & rp) const {
    return param::zero_potential_poison_value;
  }


  inline point_t
  acceleration(const body& b) const { return point_t(0); }

};

/**
  * @brief The variant captures all force types.
  *
  * @todo this could proably be automated with macros
  */

using force_var = std::variant<
    force_square_well<0>,
    force_square_well<1>,
    force_square_well<2>,
    force_spherical_wall,
    force_spherical_density_support,
    force_ndim_density_support,
    force_airfoil,
    force_gravity,
    force_central_mass,
    force_orbit,
    force_poison
  >;  // force_var


/**
  * @brief The vector of user-selectable forces
  *
  * @todo I don't love this being global, maybe the selection function should return
  *       a vector and pass ownership to the caller (likely the app driver)
  */
std::vector<force_var> vec_forces;


/**
 * @brief      Total external potential
 * @param      coords  Coordinates of where to compute the potential
 */
double
potential(const point_t & coords) {
  double phi = 0.0;
  for(auto &f : vec_forces)
    phi += std::visit([coords](auto&& x){ return x.potential(coords);}, f);
  return phi;
}

/**
 * @brief      Total external force at a point 'srch'
 * @param      particle  Accelerated particle
 */
template<class Body>
point_t
acceleration(const Body & particle) {
  point_t a = 0.0;
  for(auto &f : vec_forces)
    a += std::visit([particle](auto&& x){ return x.acceleration(particle); }, f);
  return a;
}


/**
 * @brief      External force selector
 * @param      efstr    ext. force string
 */
void
select(const std::string & efstr) {

  vec_forces.clear();

  if(boost::iequals(efstr, "zero") or boost::iequals(efstr, "none"))
    return; // trivial case

  // parse efstr: external force specification string is a comma-separated
  // list of potentials / accelerations which need to be added up: e.g.
  // "spherical wall,walls:xyz,gravity"
  std::vector<std::string> split_efstr;
  boost::split(split_efstr, efstr, boost::is_any_of(","));
  for(auto it = split_efstr.begin(); it != split_efstr.end(); ++it) {
    if(boost::iequals(*it, "spherical wall")) {
      vec_forces.emplace_back(force_spherical_wall{});
    }
    else if(boost::iequals(*it, "airfoil")) {
      vec_forces.emplace_back(force_airfoil{});
    }
    else if(boost::iequals(*it, "general density support")) {
      density_profiles::select();
      vec_forces.emplace_back(force_ndim_density_support{});
    }
    else if(boost::iequals(*it, "spherical density support")) {
      density_profiles::select();
      vec_forces.emplace_back(force_spherical_density_support{});
    }
    else if(boost::iequals(*it, "gravity")) {
      vec_forces.emplace_back(force_gravity{});
    }
    else if(boost::iequals(*it, "central mass")) {
      vec_forces.emplace_back(force_central_mass{});
    }
    else if(boost::iequals(*it, "orbit")) {
      vec_forces.emplace_back(force_orbit{});
    }
    else if(boost::iequals(it->substr(0, 6), "walls:")) {
      // parse in which directions to place the walls
      // this can be e.g. "walls:xyz" or "walls:y" etc.
      const char * ptr_cxyz = it->substr(6).c_str();
      char imx = std::min(3, (int)it->substr(6).length());
      char cxyz[3];
      strcpy(cxyz, it->substr(6, 6+imx+1).c_str());
      for(int i = 0; i < imx; ++i) {
        switch(cxyz[i]) {
          case 'x':
          case 'X':
            vec_forces.emplace_back(force_square_well<0>{});
            break;
          case 'y':
          case 'Y':
            vec_forces.emplace_back(force_square_well<1>{});
            break;
          case 'z':
          case 'Z':
            vec_forces.emplace_back(force_square_well<2>{});
            break;
          default:
            log_fatal("ERROR: bad external_force_type" << std::endl);
            assert(false);
        }
      }
    }
    else if(boost::iequals(*it, "poison")) {
      // zero potential shift
      vec_forces.emplace_back(force_poison{});
    }
    #if 0
    else if(boost::iequals(*it, "APM")) {
      vec_accelerations.push_back(artificial_pressure);
    }
    #endif
    else {
      log_fatal("ERROR: bad external_force_type" << std::endl);
    }
  } // for it in split_efstr

} // select()

/**
 * @brief      Artificial drag force - used for
 *             particle relaxation
 * @param      vel   Velocity against the drag
 */
point_t
acceleration_drag(const point_t & vel) {
  using namespace param;
  point_t acc = 0.0;
  double v2 = vel[0] * vel[0];
  for(short int i = 1; i < gdimension; ++i)
    v2 += vel[i] * vel[i];

  acc -= (relaxation_beta + relaxation_gamma * v2) * vel;
  return acc;
}


} // namespace external_force

#endif // _eforce_h_
