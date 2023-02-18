#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

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
		// shiftEvent();
		shiftEvent_efficient();

		// Return shifted results.
		allParticles_in	= allParticles;

		number_of_shifted_events = 1;

		return;
	}

	shifter::~shifter()
	{
		//clear everything

		return;
	}

	void shifter::update_records( vector<ParticleRecord> & allParticles_in )
	{
		allParticles.clear();
		pairs_sorted_by_qz.clear();
		pairs_sorted_by_abs_qz.clear();

		allParticles = allParticles_in;

		// Perform shifts.
		shiftEvent();

		// Return shifted results.
		allParticles_in	= allParticles;

		number_of_shifted_events++;

		return;
	}

	double shifter::standard_deviation( const vector<double> & v )
	{
		const double n = v.size();
		double x2 = 0.0, runsum = 0.0, variance = 0.0;
		for (const auto e: v)
		{
			variance -= e*runsum;
			runsum += e;
			x2 += e*e;
		}
		variance *= 2.0/(n*(n-1.0));
		variance += x2/n;
		return sqrt(variance);
	}



	double shifter::get_probability( const double R, const vector<double> & pair_qzs )
	{
		// if (true) return 1.0;
		/*double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		for (const auto & qz: pair_qzs) result += normalization*exp(-0.5*qz*qz*R*R);
		if (pair_qzs.size() == 3)
		{
			// extra term from Zajc's paper
			result += 2.0*exp(-0.25*R*R*(pair_qzs[0]*pair_qzs[0]
																		+pair_qzs[1]*pair_qzs[1]
																		+pair_qzs[2]*pair_qzs[2]));
		}*/
		double result = 0.0;
		double normalization = paraRdr->getVal("shifter_norm");
		for (const auto & qz: pair_qzs) result += log( 1.0 + normalization*exp(-0.5*qz*qz*R*R) );
		return exp(result);
	}

	//--------------------------------------------------------------------------
	// Perform Bose-Einstein corrections on an event.

	void shifter::shiftEvent()
	{
		//auto start = std::chrono::system_clock::now();

		// set up distributions
		bool RNG_seed	= (bool)paraRdr->getVal("RNG_seed");

		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		default_random_engine generator;
		if ( RNG_seed )
			generator = default_random_engine (seed);


	  // const int n_burn_in = 10000;
	  const double delta = 0.5;

	  std::uniform_real_distribution<double> uniform_delta(-delta, delta);
	  std::uniform_real_distribution<double> uniform(0.0, 1.0);
	  std::normal_distribution<double>       normal(0.0, 1.0);

		// burn in distribution?

		const double R = paraRdr->getVal("RNG_R") / HBARC;
		const double RNG_p0 = paraRdr->getVal("RNG_p0");
		const int number_of_particles = allParticles.size();
		vector<ParticleRecord> allParticles_Original = allParticles;

		constexpr bool check_number_of_shifted_particles = false;
		int iLoop = 0;
		// int nLoops = 1000;
		int nLoops = paraRdr->getVal("shifter_nLoops");
		for (iLoop = 0; iLoop < nLoops; iLoop++)
		{
			int number_of_shifted_particles = 0;
			for (int iParticle = 0; iParticle < number_of_particles; iParticle++)
			{
				vector<ParticleRecord> allParticles_copy = allParticles;

				// generate a shifted momentum
				// double x = allParticles.at(iParticle).p.pz();
				double x = allParticles[iParticle].p.pz();
				double y = RNG_p0 * normal(generator);	// corresponds to choice of parameters in random_events.h
				// allParticles_copy.at(iParticle).p.pz(y);
				allParticles_copy[iParticle].p.pz(y);

				// get probability of current configuration
				vector<double> unshifted_pairs = get_pairs( allParticles );
				double P1 = get_probability( R, unshifted_pairs );


				// get probability of shifted configuration
				vector<double> shifted_pairs = get_pairs( allParticles_copy );
				double P2 = get_probability( R, shifted_pairs );

				// choose new configuration (shifted or original)
		    double log_alpha = std::min(0.0, log(P2/P1) + 0.0*(y*y - x*x));
				bool shift_this_particle = (log(uniform(generator)) <= log_alpha);
		    x = shift_this_particle ? y : x;

				if (shift_this_particle) number_of_shifted_particles++;

				// store (possibly new) x to appropriate particle
				// allParticles.at(iParticle).p.pz(x);
				allParticles[iParticle].p.pz(x);
			}
			if ( check_number_of_shifted_particles && number_of_shifted_particles == 0 ) break;
		}

		// cout << "--------------------------------------------------------------------------------" << endl;
		// for (int iParticle = 0; iParticle < number_of_particles; iParticle++)
		// 	cout << "OUT: " << number_of_shifted_events << "   " << iParticle
		// 			<< "   " << allParticles_Original.at(iParticle).p.pz()
		// 			<< "   " << allParticles.at(iParticle).p.pz()
		// 			<< "   " << allParticles_Original.at(iParticle).p.pz()
		// 					- allParticles.at(iParticle).p.pz() << endl;
		vector<double> old_pairs = get_pairs( allParticles_Original );
		vector<double> new_pairs = get_pairs( allParticles );

		for (int iPair = 0; iPair < new_pairs.size(); iPair++)
			cout << "OUT: " << number_of_shifted_events << "   " << iPair
					<< "   " << old_pairs[iPair]
					<< "   " << new_pairs[iPair]
					<< "   " << new_pairs[iPair] - old_pairs[iPair] << "\n";

		cerr << "CHECK: " << standard_deviation( old_pairs )
					<< "   " << standard_deviation( new_pairs ) << "   "
					<< "   " << standard_deviation( new_pairs )
												- standard_deviation( old_pairs ) << "\n";

		cerr << "Finished shifting in " << iLoop << " of " << nLoops << " loops.\n";

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
		// Vec4 q = particles.at(i1).p - particles.at(i2).p;
		Vec4 q = particles[i1].p - particles[i2].p;
		result.push_back( q.pz() );
	}

	return result;
}


//--------------------------------------------------------------------------
// Perform Bose-Einstein corrections on an event.

void shifter::shiftEvent_efficient()
{
	// set up distributions
	bool RNG_seed	= (bool)paraRdr->getVal("RNG_seed");

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator;
	if ( RNG_seed )
		generator = default_random_engine (seed);

	std::uniform_real_distribution<double> uniform(0.0, 1.0);
	std::normal_distribution<double>       normal(0.0, 1.0);

	const double R = paraRdr->getVal("RNG_R") / HBARC;
	const double RNG_p0 = paraRdr->getVal("RNG_p0");
	const int number_of_particles = allParticles.size();
	vector<ParticleRecord> allParticles_Original = allParticles;

	vector<double> current_pairs = get_pairs(allParticles);

	constexpr bool check_number_of_shifted_particles = false;
	int iLoop = 0;
	int nLoops = paraRdr->getVal("shifter_nLoops");
	for (iLoop = 0; iLoop < nLoops; iLoop++)
	{
		int number_of_shifted_particles = 0;
		for (int iParticle = 0; iParticle < number_of_particles; iParticle++)
		{
			// generate a shifted momentum
			double x = allParticles[iParticle].p.pz();
			double y = RNG_p0 * normal(generator);	// corresponds to choice of parameters in random_events.h

			// get probability of current configuration
			double P1 = get_probability( R, current_pairs );

			// get probability of shifted configuration
			vector<ParticleRecord> allParticles_with_shift = allParticles;
			allParticles_with_shift[iParticle].p.pz(y);
			vector<double> shifted_pairs = get_shifted_pairs( current_pairs,
																			allParticles_with_shift, iParticle );
			double P2 = get_probability( R, shifted_pairs );

			// choose new configuration (shifted or original)
			double log_alpha = std::min(0.0, log(P2/P1) + 0.0*(y*y - x*x));
			bool shift_this_particle = (log(uniform(generator)) <= log_alpha);
			x = shift_this_particle ? y : x;

			if (shift_this_particle)
			{
				current_pairs = shifted_pairs;
				number_of_shifted_particles++;
			}

			// store (possibly new) x to appropriate particle
			allParticles[iParticle].p.pz(x);
		}
		if ( check_number_of_shifted_particles && number_of_shifted_particles == 0 ) break;
	}

	vector<double> old_pairs = get_pairs( allParticles_Original );
	vector<double> new_pairs = get_pairs( allParticles );

	for (int iPair = 0; iPair < new_pairs.size(); iPair++)
		cout << "OUT: " << number_of_shifted_events << "   " << iPair
				<< "   " << old_pairs[iPair]
				<< "   " << new_pairs[iPair]
				<< "   " << new_pairs[iPair] - old_pairs[iPair] << "\n";

	cerr << "CHECK: " << standard_deviation( old_pairs )
				<< "   " << standard_deviation( new_pairs ) << "   "
				<< "   " << standard_deviation( new_pairs )
											- standard_deviation( old_pairs ) << "\n";

	cerr << "Finished shifting in " << iLoop << " of " << nLoops << " loops.\n";

	// Done.
	return;

}

//--------------------------------
vector<double> shifter::get_shifted_pairs( const vector<double> & pairs,
 																					 const vector<ParticleRecord> & particles,
																					 const int shifted_particle_index )
{
	vector<double> result = pairs;
	// auto indexer = [nCols](size_t i, size_t j) { return i*nCols+j; };

	const int number_of_particles = particles.size();

	int iPair = 0;
	for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
	for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
	{
		if (i1==shifted_particle_index || i2==shifted_particle_index)
			result[iPair] = particles[i1].p.pz() - particles[i2].p.pz();
		iPair++;
	}

	return result;
}



	//--------------------------------
	bool shifter::setSortedPairs( const vector<ParticleRecord> & particles_to_sort )
	{
		// Reset.
		pairs_sorted_by_qz.clear();
		pairs_sorted_by_abs_qz.clear();

		const int number_of_particles = particles_to_sort.size();

		// get all values in vector first
		for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
		for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
		{
			// Vec4 q = particles_to_sort.at(i1).p - particles_to_sort.at(i2).p;
			Vec4 q = particles_to_sort[i1].p - particles_to_sort[i2].p;
			pairs_sorted_by_qz.push_back( std::make_pair( q.pz(), std::make_pair(i1, i2) ) );
			pairs_sorted_by_abs_qz.push_back( std::make_pair( abs(q.pz()), std::make_pair(i1, i2) ) );
		}

		// check if there are enough pairs of this species to do shift
		if (pairs_sorted_by_qz.size() < 2)
		{
			//cout << "Not enough sorted pairs of species with pid=" << idNow << endl;
			return false;
		}

		// THEN sort them (sorts on first column in ascending order automatically)
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

			// Vec4 xDiff = ( allParticles.at(i1).x - allParticles.at(i2).x ) / HBARC;
			Vec4 xDiff = ( allParticles[i1].x - allParticles[i2].x ) / HBARC;
			const double Delta_z = xDiff.pz();

			effective_source += cos(qz*Delta_z);

			npairs_in_average++;
		}

		return ( 1.0 + effective_source / npairs_in_average );
	}

}	//End of namespace

//End of file
