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
#include "math.h"

namespace external_force {

template<class Derived>
struct force_base
{

  inline double potential(const point_t& p)
  {
    return static_cast<Derived>(*this).potential(p);
  }


  inline point_t acceleration(const body& b)
  {
    return static_cast<Derived>(*this).acceleration(b);
  }
  private:
    force_base() = default;
    friend Derived;
};

template<auto I = 0>
struct force_square_well : public force_base<force_square_well<I>>
{
  const double box[3] = {.5 * param::box_length, .5 * param::box_width,
      .5 * param::box_height};


  inline double potential(const point_t& p) const {
    double phi = (((p[I] < -box[I]) ? pow(-p[I] - box[I], param::extforce_wall_powerindex) : 0.0) +
                  ((p[I] > box[I]) ? pow(p[I] - box[I], param::extforce_wall_powerindex) : 0.0)) *
                param::extforce_wall_steepness;
    return phi;
  }


  inline point_t acceleration(const body & b ) const {
    point_t a = 0.0;
    point_t rp = b.coordinates();

    a[I] = (((rp[I] < -box[I]) ? pow(-rp[I] - box[I], param::extforce_wall_powerindex - 1) : 0.0) -
            ((rp[I] > box[I]) ? pow(rp[I] - box[I], param::extforce_wall_powerindex - 1) : 0.0)) *
          param::extforce_wall_powerindex * param::extforce_wall_steepness;
    return a;
   }
};

struct force_spherical_wall : public force_base<force_spherical_wall>
{

  inline double potential(const point_t & rp) const {
    double phi = 0.0;
    double r = rp[0] * rp[0];
    for(unsigned short i = 1; i < gdimension; ++i)
      r += rp[i] * rp[i];
    r = sqrt(r);
    if(r > param::sphere_radius)
      phi = param::extforce_wall_steepness * pow(r - param::sphere_radius, param::extforce_wall_powerindex);
    return phi;
  }


  inline point_t acceleration(const body & particle) const{
    point_t a = 0.0;
    point_t rp = particle.coordinates();
    double r = rp[0] * rp[0];
    for(unsigned short i = 1; i < gdimension; ++i)
      r += rp[i] * rp[i];
    r = sqrt(r);
    if(r > param::sphere_radius) {
      const double ar = param::extforce_wall_powerindex * param::extforce_wall_steepness * pow(r - param::sphere_radius, param::extforce_wall_powerindex - 1);
      for(unsigned short i = 0; i < gdimension; ++i)
        a[i] = -rp[i] / r * ar;
    }
    return a;
  }
};

struct force_spherical_density_support : public force_base<force_spherical_density_support>
{
  const double K0 = param::pressure_initial / pow(param::rho_initial, param::poly_gamma);
  const double rho0 = density_profiles::spherical_density_profile(0.);

  force_spherical_wall _fpw;


  inline double potential(const point_t & rp) const {
    double r = rp[0] * rp[0];
    for(unsigned short i = 1; i < gdimension; ++i)
      r += rp[i] * rp[i];
    r = sqrt(r);
    const double x = r / param::sphere_radius;
    double rho =
      param::rho_initial / rho0 * density_profiles::spherical_density_profile(x);
    double phi =
      (rho > 0)
        ? (-K0 * param::poly_gamma * pow(rho, param::poly_gamma - 1.) / (param::poly_gamma - 1.))
        : 0;
    return phi + _fpw.potential(rp);
  }


  inline point_t acceleration(const body & particle) const {
    point_t a = 0.0;
    point_t rp = particle.coordinates();
    double r = rp[0] * rp[0];
    for(unsigned short i = 1; i < gdimension; ++i)
      r += rp[i] * rp[i];
    r = sqrt(r);
    const double x = r / param::sphere_radius;
    if(x > 1e-12) {
      double rho =
        param::rho_initial / rho0 * density_profiles::spherical_density_profile(x);
      double drhodr = param::rho_initial / (rho0 * param::sphere_radius) *
                      density_profiles::spherical_drho_dr(x);
      double a_r =
        (rho > 0) ? (K0 * param::poly_gamma * pow(rho, param::poly_gamma - 2) * drhodr) : 0;
      for(short int i = 0; i < gdimension; ++i)
        a[i] = a_r * rp[i] / r;
    }
    return a + _fpw.acceleration(particle);
  }

};

struct force_gravity : public force_base<force_gravity>
{

  inline double potential(const point_t & rp) const {
    double height = rp[0];
    if(gdimension > 1)
      height = rp[1];
    return height * param::gravity_acceleration_constant;;
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

struct force_airfoil : public force_base<force_airfoil>
{
  double alpha = param::airfoil_attack_angle * M_PI / 180.0;


  inline double potential(const point_t & rp) const {
    double phi = 0.0;
    assert(gdimension > 1);

    const double x1 = rp[0] - param::airfoil_anchor_x;
    const double y1 = rp[1] - param::airfoil_anchor_y;
    const double x = x1 * cos(alpha) + y1 * sin(alpha),
                y = -x1 * sin(alpha) + y1 * cos(alpha);

    bool inside_bounding_box = std::abs(y) < 5.0 * param::airfoil_thickness &&
                              x > -param::airfoil_size * 0.02 &&
                              x < param::airfoil_size * 1.02;
    double upper_surface =
      param::airfoil_thickness * x * sqrt(param::airfoil_size * param::airfoil_size - x * x);
    double camber_line = param::airfoil_camber * sin(M_PI * x / 2.);
    double aux = SQ(upper_surface) - SQ(y - camber_line) + 0.002;
    if(inside_bounding_box && aux > 0.0)
      phi = param::extforce_wall_steepness * pow(aux, param::extforce_wall_powerindex);
    return phi;
  }


  inline point_t acceleration(const body & particle) const {
    point_t a = 0.0;
    assert(gdimension > 1);

    point_t rp = particle.coordinates();
    const double x1 = rp[0] - param::airfoil_anchor_x;
    const double y1 = rp[1] - param::airfoil_anchor_y;

    const double x = x1 * cos(alpha) + y1 * sin(alpha),
                y = -x1 * sin(alpha) + y1 * cos(alpha);

    bool inside_bounding_box = std::abs(y) < 5.0 * param::airfoil_thickness &&
                              x > -param::airfoil_size * 0.02 &&
                              x < param::airfoil_size * 1.02;
    double upper_surface =
      param::airfoil_thickness * x * sqrt(param::airfoil_size * param::airfoil_size - x * x);
    double camber_line = param::airfoil_camber * sin(M_PI * x / 2.);
    double phi = SQ(upper_surface) - SQ(y - camber_line) + 0.002;
    if(inside_bounding_box && phi > 0.0) {
      double a0, a1;
      a0 = param::extforce_wall_powerindex * param::extforce_wall_steepness * pow(phi, param::extforce_wall_powerindex - 1) *
          (2. * (y - camber_line) *
              (-param::airfoil_camber * M_PI / 2. * cos(M_PI / 2. * x)) -
            param::airfoil_thickness * param::airfoil_thickness * 2 * x *
              (param::airfoil_size * param::airfoil_size - 2 * x * x));
      a1 = param::extforce_wall_powerindex * param::extforce_wall_steepness * pow(phi, param::extforce_wall_powerindex - 1) * 2. * (y - camber_line);
      a[0] = a0 * cos(alpha) - a1 * sin(alpha);
      a[1] = a0 * sin(alpha) + a1 * cos(alpha);
    }
    return a;
  }
};

struct
force_orbit : public force_base<force_orbit>
{

  double m_t = param::mass_neutron_star + param::mass_white_dwarf;


  inline double potential(const point_t & rp) const {
    assert(gdimension > 1);
    double phi = 0.0;
    // static const double grav = gravitational_constant, a_sp = orbital_separation,
    //                     m_ns = mass_neutron_star, m_wd = mass_white_dwarf;
    double term1 = sqrt(SQ(rp[0] - param::orbital_separation) + SQ(rp[1]) + SQ(rp[2]));
    term1 = -param::gravitational_constant * param::mass_neutron_star / term1;
    double term2 = -0.5 * param::gravitational_constant * m_t / CU(param::orbital_separation);
    term2 = term2 * (SQ(rp[0] - param::orbital_separation * param::mass_neutron_star / m_t) + SQ(rp[1]));
    phi = term1 + term2;
    return phi;
  }


  inline point_t acceleration(const body & particle) const {
    point_t rp = particle.coordinates();
    point_t acc = 0.0;

    double temp = SQ(rp[0] - param::orbital_separation) + SQ(rp[1]) + SQ(rp[2]);
    temp = CU(temp);
    temp = sqrt(temp);
    double term1 = -param::gravitational_constant * param::mass_neutron_star / temp;
    acc[0] += term1 * (rp[0] - param::orbital_separation);
    acc[1] += term1 * rp[1];
    acc[2] += term1 * rp[2];

    double term2 = param::gravitational_constant * m_t / CU(param::orbital_separation);
    acc[0] += term2 * (rp[0] - param::orbital_separation * param::mass_neutron_star / m_t); // x-direction
    acc[1] += term2 * rp[1];
    return acc;
  }

};

struct force_poison : public force_base<force_poison>
{

  inline double
  potential(const point_t & rp) const {
    return param::zero_potential_poison_value;
  }


  inline point_t
  acceleration(const body& b) const { return point_t(0); }

};

using force_var = std::variant< force_square_well<0>, force_square_well<1>, force_square_well<2>,
                                force_spherical_wall, force_spherical_density_support, force_airfoil,
                                force_gravity, force_orbit, force_poison>;

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
    else if(boost::iequals(*it, "spherical density support")) {
      density_profiles::select();
      vec_forces.emplace_back(force_spherical_density_support{});
    }
    else if(boost::iequals(*it, "gravity")) {
      vec_forces.emplace_back(force_gravity{});
    }
    else if(boost::iequals(*it, "orbit")) {
      vec_forces.emplace_back(force_orbit{});
    }
    else if(boost::iequals(it->substr(0, 6), "walls:")) {
      // parse in which directions to place the walls
      // this can be e.g. "walls:xyz" or "walls:y" etc.
      const char * cxyz = it->substr(6).c_str();
      char imx = std::min(3, (int)it->substr(6).length());
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
