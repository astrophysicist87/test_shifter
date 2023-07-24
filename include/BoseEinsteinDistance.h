#ifndef BOSEEINSTEINDISTANCE_H
#define BOSEEINSTEINDISTANCE_H

#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

#include "param_list.h"
#include "ParameterReader.h"
#include "Particle.h"

using shift_lib::Particle;
using shift_lib::particle_pair;

//==============================================================================
class BoseEinsteinDistance
{
  //----------------------------------------------------------------------------
  private:
    double R{0.0}, m{0.13957};
		ostream & out;
		ostream & err;

    inline double norm2(std::initializer_list<double> il)
    {
      return inner_product( il.begin(), il.end(), il.begin(), 0.0 );
    }

    inline double get_distance_single_scale(const Particle & p1, const Particle & p2)
  	{
      // assume units are fine
      double adhoc_R2 = R*R/(1.0
                            + 3.0*sqrt(m*m
                                      + 0.25*norm2( { p1.p.x() + p2.p.x(),
                                                      p1.p.y() + p2.p.y(),
                                                      p1.p.z() + p2.p.z() }) ));
      // note: 0.25 should eventually be replaced by 0.5
      return exp( -0.25 * adhoc_R2 * norm2( { p1.p.x() - p2.p.x(),
                                              p1.p.y() - p2.p.y(),
                                              p1.p.z() - p2.p.z() } ) );
      // return exp( -0.25 * R * R * norm2( { p1.p.x() - p2.p.x(),
      //                                      p1.p.y() - p2.p.y(),
      //                                      p1.p.z() - p2.p.z() } ) );
      // return exp( -0.25 * R * R * norm2( { p1.p.x() - p2.p.x(),
      //                                      p1.p.y() - p2.p.y(),
      //                                      p1.p.z() - p2.p.z() } )
      //             -0.125 * norm2( { p1.p.x() + p2.p.x(),
      //                                      p1.p.y() + p2.p.y(),
      //                                      p1.p.z() + p2.p.z() } ) );
    }

  //----------------------------------------------------------------------------
  public:
    BoseEinsteinDistance( const string & mode,
                  				const param_list & parameters,
                  				ostream & out_stream = std::cout,
                  				ostream & err_stream = std::cerr )
		: out{out_stream},
		  err{err_stream}
		{
			if (mode == "SingleScale")
			{
				R = std::get<double>(parameters.at("R"));
        get_distance = [this](const Particle & p1, const Particle & p2)
                       { return get_distance_single_scale(p1, p2); };
			}
			else
			{
				err << "BoseEinsteinDistance::Invalid mode = " << mode << endl;
				std::terminate();
			}
		}
    ~BoseEinsteinDistance(){}

    std::function<double(const Particle &, const Particle &)> get_distance;
};

#endif
