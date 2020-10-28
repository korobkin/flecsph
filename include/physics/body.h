/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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
 * @file body.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Representation of a body, a particle for our SPH implementation
 */

#ifndef body_h
#define body_h

#define OUTPUT

#include "space_vector.h"
#include "tree_topology/tree_types.h"
#include "phys_consts.h"
#include "user.h"

enum particle_type_t : int { NORMAL = 0, WALL = 1 };

enum state_t : int { NONE = 0, STAR1 = 1, STAR2 = 2, POINTP = 3 };

template<class KEY>
class body_u : public flecsi::topology::entity<gdimension, type_t, KEY>
{

  static const size_t dimension = gdimension;
  using element_t = type_t;
  using point_t = flecsi::space_vector_u<element_t, dimension>;

  using flecsi::topology::entity<gdimension, type_t, KEY>::mass_;

public:
  body_u()
    : flecsi::topology::entity<gdimension, type_t, KEY>(), type_(NORMAL),
      state_(NONE){};

  double getPressure() const {
    return pressure_;
  }
  double getSoundspeed() const {
    return soundspeed_;
  }
  double getEntropy() const {
    return entropy_;
  }
  double getElectronfraction() const {
    return electronfraction_;
  }
  double getAbar() const {
    return abar_;
  }
  double getDensity() const {
    return density_;
  }
  double getTemperature() const {
    return temperature_;
  }
  point_t getVelocity() const {
    return velocity_;
  }
  point_t getVelocityhalf() const {
    return velocityhalf_;
  }
  point_t getAcceleration() const {
    return acceleration_;
  }
  point_t getGAcceleration() const {
    return g_acceleration_;
  }
  particle_type_t type() const {
    return type_;
  };
  double getGPotential() const {
    return g_potential_;
  }

  point_t getLinMomentum() const {
    point_t res = {};
    for(size_t i = 0; i < dimension; ++i) {
      res[i] = velocity_[i] * mass_;
    }
    return res;
  };
  double getDt() const {
    return dt_;
  };
  particle_type_t getType() const {
    return type_;
  };
  state_t state() const {
    return state_;
  }
  bool is_wall() const {
    return type_ == 1;
  };

  // Accessors "get<Function>InGeom() return corresponding 
  // quantities converted to geometric units
  double getPressureInGeom() const {
    return 1./P_GEOM_TO_CGS*pressure_;
  }
  double getDensityInGeom() const {
    return 1./RHO_GEOM_TO_CGS*density_;
  }
  double getInternalenergyInGeom() const {
    return 1./EPS_GEOM_TO_CGS*internalenergy_;
  }
  double getSoundspeedInGeom() const {
    return 1./VEL_GEOM_TO_CGS*soundspeed_;
  }
  double getMassInGeom() const {
    return 1./M_GEOM_TO_CGS*this->mass();
  }
  double getRadiusInGeom() const {
    return 1./L_GEOM_TO_CGS*this->radius();
  }
  point_t getCoordinatesInGeom() const {
    return 1./L_GEOM_TO_CGS*this->coordinates();
  }
  point_t getVelocityInGeom() const {
    return 1./VEL_GEOM_TO_CGS*velocity_;
  }
  point_t getVelocityhalfInGeom() const {
    return 1./VEL_GEOM_TO_CGS*velocityhalf_;
  }
  point_t getAccelerationInGeom() const {
    return 1./ACC_GEOM_TO_CGS*acceleration_;
  }
 
  // --- Mutators "set<Function>InGeom":
  //     take an argument in geometric units,
  //     convert it to CGS and write the result in CGS,
  //     so ultimately everything is stored in CGS
  void setPressureInGeom(const double pressure) {
    pressure_ = P_GEOM_TO_CGS*pressure;
  }
  void setDensityInGeom(const double density) {
    density_ = RHO_GEOM_TO_CGS*density;
  }
  void setInternalenergyInGeom(const double internalenergy) {
    internalenergy_ = EPS_GEOM_TO_CGS*internalenergy;
  }
  void setSoundspeedInGeom(const double soundspeed) {
    soundspeed_ = VEL_GEOM_TO_CGS*soundspeed;
  }
  void setMassInGeom(const double mass) {
    this->set_mass(M_GEOM_TO_CGS*mass);
  }
  void setRadiusInGeom(const double radius) {
    this->set_radius(L_GEOM_TO_CGS*radius);
  }
  void setCoordinatesInGeom(const point_t & coordinates) {
    this->set_coordinates(L_GEOM_TO_CGS*coordinates);
  }
  void setVelocityInGeom(const point_t & velocity) {
    velocity_ = VEL_GEOM_TO_CGS*velocity;
  }
  void setVelocityhalfInGeom(const point_t & velocityhalf) {
    velocityhalf_ = VEL_GEOM_TO_CGS*velocityhalf;
  }
  void setAccelerationInGeom(const point_t & acceleration) {
    acceleration_ = ACC_GEOM_TO_CGS*acceleration;
  }
  // ---

  void setAcceleration(const point_t & acceleration) {
    acceleration_ = acceleration;
  }
  void setGAcceleration(const point_t & g_acceleration) {
    g_acceleration_ = g_acceleration;
  }
  void setGPotential(const double g_potential) {
    g_potential_ = g_potential;
  }
  void setVelocity(const point_t & velocity) {
    velocity_ = velocity;
  }
  void setVelocityhalf(const point_t & velocityhalf) {
    velocityhalf_ = velocityhalf;
  }
  void setSoundspeed(const double soundspeed) {
    soundspeed_ = soundspeed;
  }
  void setPressure(const double pressure) {
    pressure_ = pressure;
  }
  void setEntropy(const double entropy) {
    entropy_ = entropy;
  }
  void setElectronfraction(const double electronfraction) {
    electronfraction_ = electronfraction;
  }
  void setAbar(const double abar) {
    abar_ = abar;
  }
  void setDensity(const double density) {
    density_ = density;
  }
  void setTemperature(const double temperature) {
    temperature_ = temperature;
  }
  void setDt(const double dt) {
    dt_ = dt;
  }
  void setType(const particle_type_t type) {
    type_ = type;
  }
  void setType(const int type) {
    type_ = static_cast<particle_type_t>(type);
  }
  void set_state(const state_t state) {
    state_ = state;
  }
  void set_state(const int state) {
    state_ = static_cast<state_t>(state);
  }
  double getInternalenergy() const {
    return internalenergy_;
  }
  void setInternalenergy(const double internalenergy) {
    internalenergy_ = internalenergy;
  }
  double getTotalenergy() const {
    return totalenergy_;
  }
  void setTotalenergy(const double totalenergy) {
    totalenergy_ = totalenergy;
  }
  void setDudt(const double dudt) {
    dudt_ = dudt;
  }
  double getDudt() const {
    return dudt_;
  }
  void setDedt(const double dedt) {
    dedt_ = dedt;
  }
  double getDedt() const {
    return dedt_;
  }
  void setAlpha(const double alpha) {
    alpha_ = alpha;
  }
  double getAlpha() const {
    return alpha_;
  }
  void setDivergenceV(const double divergenceV) {
    divergenceV_ = divergenceV;
  }
  double getDivergenceV() const {
    return divergenceV_;
  }
  void setDdivvdt(const double dDivVdt) {
    dDivVdt_ = dDivVdt;
  }
  double getDdivvdt() const {
    return dDivVdt_;
  }
  void setTrigger(const double trigger) {
    trigger_ = trigger;
  }
  double getTrigger() const {
    return trigger_;
  }
  void setXi(const double xi) {
    xi_ = xi;
  }
  double getXi() const {
    return xi_;
  }
  void setTraceSS(const double traceSS) {
    traceSS_ = traceSS;
  }
  double getTraceSS() const {
    return traceSS_;
  }

  void setGradV(const double gradv) {
    gradv_ = gradv;
  }
  double getGradV() const {
    return gradv_;
  }

  void setNeighbors(const size_t neighbors) {
    neighbors_ = neighbors;
  }
  size_t getNeighbors() const {
    return neighbors_;
  }

  void setSignalspeed(const double & signalspeed) {
    signalspeed_ = signalspeed;
  }
  double getSignalspeed() const {
    return signalspeed_;
  }

  friend std::ostream & operator<<(std::ostream & os, const body_u & b) {
    // TODO change regarding to dimension
    os << std::setprecision(10);
    os << "Particle: Pos: " << b.coordinates_ << " rho: " << b.density_;
    os << " h: " << b.radius_;
    os << " P: " << b.pressure_;
    os << " v: " << b.velocity_;
    os << " m: " << b.mass_;
    os << " u: " << b.internalenergy_;
    os << " cs: " << b.soundspeed_;
    os << " a: " << b.acceleration_;
    os << " ga: " << b.g_acceleration_;
    os << "gpot: " << b.g_potential_;
    os << " id: " << b.id_;
    os << " key: " << b.key_;
    os << " owner: " << b.owner_;
    os << " neighbors: " << b.neighbors_;
    return os;
  }

private:
  point_t velocity_;
  point_t velocityhalf_;
  point_t acceleration_;
  point_t g_acceleration_;
  double g_potential_;
  double density_;
  double pressure_;
  double entropy_;
  double electronfraction_;
  double abar_;
  double temperature_;
  double soundspeed_;
  double internalenergy_;
  double totalenergy_;
  double dudt_;
  double dedt_;
  double dt_;
  double alpha_;
  double divergenceV_;
  double dDivVdt_;
  double trigger_;
  double xi_;
  double traceSS_;
  double gradv_;
  particle_type_t type_;
  size_t neighbors_;
  state_t state_;
  double signalspeed_;
}; // class body

#endif // body_h
