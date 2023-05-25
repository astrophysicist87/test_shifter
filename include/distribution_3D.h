#ifndef DISTRIBUTION_3D_H
#define DISTRIBUTION_3D_H

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <variant>

#include "param_list.h"

using namespace std;

class distribution_3D
{

private:

	bool do_x = false, do_y = false, do_z = false;
	double mean = 0.0, scale = 0.0, width = 0.0;

	default_random_engine generator;
	normal_distribution<double> normal_dist;
	uniform_real_distribution<double> pMag_dist;
	uniform_real_distribution<double> pTheta_dist;
	uniform_real_distribution<double> pPhi_dist;

	ostream & out;
	ostream & err;

	std::tuple<double,double,double> generate_sample_Gaussian()
	{
		double px = do_x ? normal_dist(generator) : 0.0;
		double py = do_y ? normal_dist(generator) : 0.0;
		double pz = do_z ? normal_dist(generator) : 0.0;
		return std::make_tuple(px,py,pz);
	}

	std::tuple<double,double,double> generate_sample_Shell()
	{
		double pMag   = pMag_dist(generator);
		double pTheta = pTheta_dist(generator);
		double pPhi   = pPhi_dist(generator);
		double sin_pTheta = sin(pTheta);
		double px = pMag*sin_pTheta*cos(pPhi);
		double py = pMag*sin_pTheta*sin(pPhi);
		double pz = pMag*cos(pTheta);
		return std::make_tuple(px,py,pz);
	}

public:

	distribution_3D( const string & mode, int seed,
				param_list & parameters,
				ostream & out_stream = std::cout,
				ostream & err_stream = std::cerr )
				:
				out(out_stream),
				err(err_stream),
				generator( default_random_engine(seed) )
				{
					if (mode == "Gaussian")
					{
						mean  = std::get<double>(parameters.at("mean"));
						scale = std::get<double>(parameters.at("scale"));
						do_x  = std::get<bool>(parameters.at("do_x"));
						do_y  = std::get<bool>(parameters.at("do_y"));
						do_z  = std::get<bool>(parameters.at("do_z"));
						normal_dist = normal_distribution<double>(mean, scale);
						generate_sample = [this]{ return generate_sample_Gaussian(); };
					}
					else if (mode == "Shell")
					{
						scale = std::get<double>(parameters.at("scale"));
						width = std::get<double>(parameters.at("width"));
						pMag_dist = uniform_real_distribution<double>(scale-width, scale+width);
						pTheta_dist = uniform_real_distribution<double>(0.0, M_PI);
						pPhi_dist = uniform_real_distribution<double>(-M_PI, M_PI);
						generate_sample = [this]{ return generate_sample_Shell(); };
					}
					else
					{
						err << "Invalid mode = " << mode << endl;
						std::terminate();
					}
				}

	~distribution_3D(){}

  std::function<std::tuple<double,double,double>()> generate_sample;

};

#endif
