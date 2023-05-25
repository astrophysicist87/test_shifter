#ifndef QUANTUMSAMPLER_H
#define QUANTUMSAMPLER_H

#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <unordered_set>
#include <vector>

#include "ParameterReader.h"
#include "Particle.h"
#include "stopwatch.h"

typedef std::uniform_int_distribution<std::mt19937::result_type> integer_range;
typedef std::normal_distribution<double>                         gaussian;

using shift_lib::Particle;

constexpr double TINY = 1e-6;
constexpr double HBARC = 0.19733;
//==============================================================================
class QuantumSampler
{
  //----------------------------------------------------------------------------
  private:
    //--------------------------------------------------------------------------
    int RNG_xDir { 0 };
    int RNG_yDir { 0 };
    int RNG_zDir { 0 };
    int dimension { 0 };
    double sigma { 1.0 };  // fm
    vector<Particle> & particles;
    std::mt19937 rng;
    default_random_engine generator;
    integer_range intdist;
    gaussian normdist;

    //--------------------------------------------------------------------------
    inline int choose_particle() { return intdist(rng); };

  //----------------------------------------------------------------------------
  public:
    QuantumSampler( vector<Particle> & particles_in,
                    double sigma_in, int RNG_xDir_in, int RNG_yDir_in, int RNG_zDir_in,
                    unsigned seed = chrono::system_clock::now().time_since_epoch().count() )
    : particles{ particles_in },
      rng{ std::mt19937(seed) },
      intdist{ integer_range(0, particles_in.size() - 1) },
      normdist{ gaussian(0.0, 1.0) }
    {
      // paraRdr = paraRdr_in;
      // sigma = paraRdr->getVal("sigma")/HBARC;  // convert to inverse GeV
      // RNG_xDir = paraRdr->getVal("RNG_xDir");
      // RNG_yDir = paraRdr->getVal("RNG_yDir");
      // RNG_zDir = paraRdr->getVal("RNG_zDir");
      sigma = sigma_in/HBARC;  // convert to inverse GeV
      RNG_xDir = RNG_xDir_in;
      RNG_yDir = RNG_yDir_in;
      RNG_zDir = RNG_zDir_in;
      dimension = RNG_xDir + RNG_yDir + RNG_zDir;
      // s1p_normalization = pow(2.0*sqrt(M_PI)*sigma, dimension);
    }
    ~QuantumSampler(){}

    std::tuple<double,double,double> sample(
      std::tuple<double,double,double> def = std::make_tuple(0.0, 0.0, 0.0))
    {
      auto & pi = particles[ choose_particle() ].p;
      double w = 1.0/(sqrt(2.0)*sigma);
      double rx = RNG_xDir ? normdist(rng) : get<0>(def);
      double ry = RNG_yDir ? normdist(rng) : get<1>(def);
      double rz = RNG_zDir ? normdist(rng) : get<2>(def);
      return std::make_tuple( pi.x()+w*rx, pi.y()+w*ry, pi.z()+w*rz );
    }
};

#endif
