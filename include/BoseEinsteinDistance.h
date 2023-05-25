#ifndef BOSEEINSTEINDISTANCE_H
#define BOSEEINSTEINDISTANCE_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;

#include "param_list.h"
#include "ParameterReader.h"

//==============================================================================
class BoseEinsteinDistance
{
  //----------------------------------------------------------------------------
  private:
    double R{0.0};
		ostream & out;
		ostream & err;

    inline double get_distance_single_scale(const vector<double> & q)
  	{
      return exp(-0.25 * R * R * (q[0]*q[0] + q[1]*q[1] + q[2]*q[2]));
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
				get_distance = [this](const vector<double> & q)
                       { return get_distance_single_scale(q); };
			}
			else
			{
				err << "BoseEinsteinDistance::Invalid mode = " << mode << endl;
				std::terminate();
			}
		}
    ~BoseEinsteinDistance(){}

    std::function<double(const vector<double> &)> get_distance;
};

#endif
