#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

#include "./ParameterReader.h"

using namespace std;

namespace shift_lib
{

	std::tuple<double,double,double> generate_sample_from_spectra_Gaussian(
		default_random_engine & generator,
		normal_distribution<double> & distribution,
		ParameterReader * paraRdr )
	{
		double RNG_p0 = paraRdr->getVal("RNG_p0");
		int RNG_xDir 	= paraRdr->getVal("RNG_xDir");
		int RNG_yDir 	= paraRdr->getVal("RNG_yDir");
		int RNG_zDir 	= paraRdr->getVal("RNG_zDir");
		double px = RNG_xDir ? RNG_p0 * distribution(generator) : 0.0;
		double py = RNG_yDir ? RNG_p0 * distribution(generator) : 0.0;
		double pz = RNG_zDir ? RNG_p0 * distribution(generator) : 0.0;
		return std::make_tuple(px,py,pz);
	}

	std::tuple<double,double,double> generate_sample_from_spectra_shell(
		default_random_engine & generator,
		normal_distribution<double> & distribution,
		ParameterReader * paraRdr )
	{
		double RNG_p0 = paraRdr->getVal("RNG_p0");
		int RNG_xDir 	= paraRdr->getVal("RNG_xDir");
		int RNG_yDir 	= paraRdr->getVal("RNG_yDir");
		int RNG_zDir 	= paraRdr->getVal("RNG_zDir");
		double px = RNG_xDir ? RNG_p0 * distribution(generator) : 0.0;
		double py = RNG_yDir ? RNG_p0 * distribution(generator) : 0.0;
		double pz = RNG_zDir ? RNG_p0 * distribution(generator) : 0.0;
		return std::make_tuple(px,py,pz);
	}
}

#endif
