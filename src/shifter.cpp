#include <iostream>
#include <fstream>
#include <ios>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>
#include <random>
#include <algorithm>

#include "../include/shifter.h"
#include "../include/ParameterReader.h"

using namespace std;

namespace shift_lib
{

	constexpr bool BE_VERBOSE = false;

	constexpr double MM2FM = 1e12;
	constexpr double FM2MM = 1e-12;
	constexpr double HBARC = 0.19733;

	constexpr double COMPRELERR = 1e-10;
	constexpr double COMPFACMAX = 1000.;
	constexpr int    NCOMPSTEP  = 10;


	void shifter::initialize_all( ParameterReader * paraRdr_in,
		vector<ParticleRecord> & allParticles_in )
	{
		include_pair_density = false;

		// Load parameters.
		paraRdr			= paraRdr_in;

		// Load particles.
		allParticles	= allParticles_in;

		// Perform shifts.
		shiftEvent();

		// Return shifted results.
		allParticles_in	= allParticles;

		return;
	}

	shifter::~shifter()
	{
		//clear everything

		return;
	}

	void shifter::update_records( vector<ParticleRecord> & allParticles_in )
	{

		return;
	}

	double shifter::get_probability( const double R, const vector<double> & pair_qzs )
	{
		double result = 1.0;
		for (const auto & qz: pairs) result += exp(-qz*qz*R*R);
		return result;
	}

	//--------------------------------------------------------------------------
	// Perform Bose-Einstein corrections on an event.

	void shifter::shiftEvent()
	{
		//auto start = std::chrono::system_clock::now();

		// set up distributions
	  std::default_random_engine generator;
	  const int n_burn_in = 10000;
	  const double delta = 0.5;

	  std::uniform_real_distribution<double> uniform_delta(-delta, delta);
	  std::normal_distribution<double> normal_delta(0.0, 1.0);
	  std::uniform_real_distribution<double> uniform(0.0, 1.0);

		// burn in distribution?

		const double R = 5.0 / HBARC;
		const int number_of_particles = allParticles.size();
		vector<double> allParticles_copy = allParticles;

		for (int iParticle = 0; iParticle < number_of_particles; iParticle++)
		{
			// generate a shifted momentum
			double x = allParticles.at(iParticle).p.pz();
    	double y = x + uniform_delta(generator);
			allParticles_copy.at(iParticle).p.pz(y);

			// get probability of current configuration
			vector<double> unshifted_pairs = get_pairs( allParticles );
			double P1 = get_probability( R, unshifted_pairs );


			// get probability of shifted configuration
			vector<double> shifted_pairs = get_pairs( allParticles_copy );
			double P2 = get_probability( R, shifted_pairs );

			// choose new configuration (shifted or original)
		}

		// Done.
		return;

	}

	//--------------------------------
	vector<double> shifter::get_pairs( const vector<ParticleRecord> & particles )
	{
		vector<double> result;

		const int number_of_particles = particles.size();

		// get all values in vector first
		for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
		for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
		{
			Vec4 q = particles.at(i1).p - particles.at(i2).p;
			result.push_back( q.pz() );
		}

		return result;
	}

	//--------------------------------
	bool shifter::setSortedPairs( const vector<ParticleRecord> & particles_to_sort )
	{
		// Reset.
		sortedPairs.clear();
		pairs_sorted_by_qz.clear();
		pairs_sorted_by_abs_qz.clear();

		const int number_of_particles = particles_to_sort.size();

		// get all values in vector first
		for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
		for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
		{
			Vec4 q = particles_to_sort.at(i1).p - particles_to_sort.at(i2).p;
			sortedPairs.push_back( std::make_pair( qPRF.pAbs(), std::make_pair(i1, i2) ) );
			pairs_sorted_by_qz.push_back( std::make_pair( q.pz(), std::make_pair(i1, i2) ) );
			pairs_sorted_by_abs_qz.push_back( std::make_pair( abs(q.pz()), std::make_pair(i1, i2) ) );
		}

		// check if there are enough pairs of this species to do shift
		if (sortedPairs.size() < 2)
		{
			//cout << "Not enough sorted pairs of species with pid=" << idNow << endl;
			return false;
		}

		sortedPairs.push_back( std::make_pair( 0.0, std::make_pair(-1, -1) ) );

		// THEN sort them (sorts on first column in ascending order automatically)
		sort( sortedPairs.begin(), sortedPairs.end() );
		sort( pairs_sorted_by_qz.begin(), pairs_sorted_by_qz.end() );
		sort( pairs_sorted_by_abs_qz.begin(), pairs_sorted_by_abs_qz.end() );

		return (true);
	}



	// current effective source is 1 + < cos(qz Delta_z) >
	double shifter::evaluate_effective_source(
					const vector< pair< double, pair <int,int> > > & pairs,
					const double qz )
	{
		int npairs_in_average = 0;
		double effective_source = 0.0;

		// the BE enhancement piece
		for (const auto & iPair : pairs)
		{
			const int i1 = iPair.second.first;
			const int i2 = iPair.second.second;

			if ( i1<0 or i2<0 ) continue;
			//if ( this1 != i1 and this2 != i2 ) continue;
			//if ( this1 != i1 or  this2 != i2 ) continue;

			Vec4 xDiff = ( allParticles.at(i1).x - allParticles.at(i2).x ) / HBARC;
			const double Delta_z = xDiff.pz();

			effective_source += cos(qz*Delta_z);

			npairs_in_average++;
		}

		return ( 1.0 + effective_source / npairs_in_average );
	}

}	//End of namespace

//End of file
