#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../include/shifter.h"
#include "../include/ParameterReader.h"

#include "KolmogorovSmirnovDist.c"

using namespace std;

namespace shift_lib
{
	constexpr const char * SHIFT_MODE = "TRIAL3";

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
		// shiftEvent();
		shiftEvent_efficient();

		// Return shifted results.
		allParticles_in	= allParticles;

		number_of_shifted_events++;

		return;
	}

	void shifter::get_combinations(int N, int K, vector<vector<int>> & combinations)
	{
		std::string bitmask(K, 1); // K leading 1's
		bitmask.resize(N, 0); // N-K trailing 0's

		// store integers and permute bitmask
		do
		{
			int ind = 0;
			vector<int> combination(3);

	    for (int i = 0; i < N; ++i) // [0..N-1] integers
        if (bitmask[i]) combination[ind++] = i;

			combinations.push_back( combination );
		} while (std::prev_permutation(bitmask.begin(), bitmask.end()));

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



	// double shifter::get_probability( const double R, const vector<double> & pair_qzs )
	// {
	// 	// if (true) return 1.0;
	// 	// double result = 1.0;
	// 	// double normalization = paraRdr->getVal("shifter_norm");
	// 	// for (const auto & qz: pair_qzs) result += normalization*exp(-0.5*qz*qz*R*R);
	// 	// if (pair_qzs.size() == 3)
	// 	// {
	// 	// 	// extra term from Zajc's paper
	// 	// 	result += 2.0*exp(-0.25*R*R*(pair_qzs[0]*pair_qzs[0]
	// 	// 																+pair_qzs[1]*pair_qzs[1]
	// 	// 																+pair_qzs[2]*pair_qzs[2]));
	// 	// }
	// 	// return result;
	// 	//--------------------------------------------------------------------------
	// 	// double result = 0.0;
	// 	// double normalization = paraRdr->getVal("shifter_norm");
	// 	// for (const auto & qz: pair_qzs) result += log( 1.0 + normalization*exp(-0.5*qz*qz*R*R) );
	// 	// return exp(result);
	// 	//--------------------------------------------------------------------------
	// 	// FULL PRODUCT
	// 	// double result = 1.0;
	// 	// double normalization = paraRdr->getVal("shifter_norm");
	// 	// for (const auto & qz: pair_qzs) result *= 1.0 + normalization*exp(-0.5*qz*qz*R*R);
	// 	// return result;
	// 	//--------------------------------------------------------------------------
	// 	// PRODUCT TO SECOND ORDER
	// 	// const int n = pair_qzs.size();
	// 	// double result = 1.0;
	// 	// double norm = paraRdr->getVal("shifter_norm");
	// 	// for (const auto & qz: pair_qzs) result += exp(-0.5*qz*qz*R*R);
	// 	// for (int ij = 0;    ij < n; ij++)
	// 	// for (int kl = ij+1; kl < n; kl++)
	// 	// 	result += norm*norm*exp(-0.5*(pair_qzs[ij]*pair_qzs[ij]+pair_qzs[kl]*pair_qzs[kl])*R*R);
	// 	// return result;
	// 	//--------------------------------------------------------------------------
	// 	// CRAMER VON MISES STATISTIC
	// 	const int n = pair_qzs.size();
	// 	double T = 1.0/(12.0*n);
	// 	for (int i = 0; i < n; i++)
	// 	{
	// 		double xi = (2.0*i+1.0)/(2.0*n);
	// 		double Fi = 0.5*(1.0+erf(pair_qzs[i]*R/sqrt(2.0)));
	// 		T += (xi-Fi)*(xi-Fi);
	// 	}
	// 	return 1.0/sqrt(T);
	// }


//------------------------------------------------------------------------------
// THIS IS ONLY FOR COMPARING PURPOSES
double shifter::get_probability( const double R, const vector<double> & pair_qzs )
{
	//----------------------------------------------------------------------------
	// KOLMOGOROV-SMIRNOV STATISTIC
	// const int n = pair_qzs.size();
	// double KS_D = -1.0;
	// for (int i = 1; i <= n; i++)
	// {
	// 	double Fi = 0.5*(1.0+erf(pair_qzs[i-1]*R/sqrt(2.0)));
	// 	KS_D = std::max( KS_D, std::max( Fi - (i-1.0)/n, (i*1.0/n) - Fi ) );
	// }
	// if (KS_D < 0.0 || KS_D > 1.0)
	// {
	// 	cerr << "KS_D = " << KS_D << endl;
	// 	abort();
	// }
	// //auto KS_n1 = [](double x) { return ( (x < 0.5) ? 0.0 : 2.0*x-1.0 ); };
	// //return 2.0-KS_n1(KS_D);
	// return 2.0-KScdf( n, KS_D );
	//----------------------------------------------------------------------------
	// MODULATE ONLY LARGEST DEVIATION FROM CDF
	// const int n = pair_qzs.size();
	// double max_deviation = -1.0;
	// int max_i = -1;
	// for (int i = 1; i <= n; i++)
	// {
	// 	double Fi = 0.5*(1.0+erf(pair_qzs[i-1]*R/sqrt(2.0)));
	// 	double this_deviation = std::max( Fi - (i-1.0)/n, (i*1.0/n) - Fi );
	// 	if ( this_deviation > max_deviation )
	// 	{
	// 		max_i = i;
	// 		max_deviation = this_deviation;
	// 	}
	// }
	// if (max_deviation < 0.0 || max_deviation > 1.0 || max_i < 1 || max_i > n)
	// {
	// 	cerr << "max_deviation = " << max_deviation << endl;
	// 	cerr << "max_i = " << max_i << endl;
	// 	abort();
	// }
	// double qz = pair_qzs[max_i-1];
	// return 1.0+exp(-0.5*qz*qz*R*R);
	//--------------------------------------------------------------------------
	// CRAMER VON MISES STATISTIC
	if ( SHIFT_MODE == "CramerVonMises" )
	{
		const int n = pair_qzs.size();
		double T = 1.0/(12.0*n);
		for (int i = 0; i < n; i++)
		{
			double xi = (2.0*i+1.0)/(2.0*n);
			double Fi = 0.5*(1.0+erf(pair_qzs[i]*R/sqrt(2.0)));
			T += (xi-Fi)*(xi-Fi);
		}
		return 1.0/sqrt(T);
	}
	//--------------------------------------------------------------------------
	// FULL PRODUCT
	else if ( SHIFT_MODE == "FullProduct" )
	{
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		for (const auto & qz: pair_qzs) result *= 1.0 + normalization*exp(-0.5*qz*qz*R*R);
		return result;
	}
	//--------------------------------------------------------------------------
	// FULL PRODUCT RAISED TO A POWER
	else if ( SHIFT_MODE == "FullProductPower" )
	{
		const int n = pair_qzs.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		for (const auto & qz: pair_qzs) result *= 1.0 + 0.5*np*normalization*exp(-0.5*qz*qz*R*R);
		return std::pow(result, 2.0/np);
	}
	//--------------------------------------------------------------------------
	// EXPERIMENTAL
	else if ( SHIFT_MODE == "TRIAL" )
	{
		const int n = pair_qzs.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		for (int i = 0; i < np-1; i++)
			result *= 1.0 + normalization*exp(-0.5*pair_qzs[i]*pair_qzs[i]*R*R);
		return result;
	}
	else if ( SHIFT_MODE == "TRIAL2" )
	{
		const int n = pair_qzs.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		int i = -1;
		for (int i1 = 0; i1 < np - 1; ++i1)
		for (int i2 = i1 + 1; i2 < np; ++i2)
		{
			i++;
			bool include_this_pair = (i2 == i1+1) /*|| (i1 == 0 && i2 == np-1)*/;
			if (!include_this_pair) continue;
			result *= 1.0 + (0.5*np/**(np-1.)/np*/)*normalization*exp(-0.5*pair_qzs[i]*pair_qzs[i]*R*R);
		}
		return result;
	}
	else if ( SHIFT_MODE == "TRIAL3" )
	{
		const int n = pair_qzs.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		int i = -1;
		double factor = 0.0;
		for (int i1 = 0; i1 < np - 1; ++i1)
		for (int i2 = i1 + 1; i2 < np; ++i2)
		{
			i++;
			bool include_this_pair = (i2 == i1+1) || (i1 == 0 && i2 == np-1);
			if (!include_this_pair) continue;
			double term = 1.0 + 0.5*np*normalization*exp(-0.5*pair_qzs[i]*pair_qzs[i]*R*R);
			result *= term;
			factor += 1.0/term;
		}
		return factor*result/np;
	}
	else if ( SHIFT_MODE == "TRIAL4" )
	{
		const int n = pair_qzs.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");

		// compute all n-choose-(np-1) combinations of pairs
		vector<vector<int>> combinations;
		get_combinations(n, np-1, combinations);

		// store all factors for each pair
		vector<double> factors(n);
		for (int i = 0; i < n; i++)
			factors[i] = 1.0 + 0.5*np*normalization*exp(-0.5*pair_qzs[i]*pair_qzs[i]*R*R);

		// sum all products of terms and return average
		double sum = 0.0;
		for ( const auto & c : combinations )
		{
			double term = 1.0;
			for ( const int & i : c ) term *= factors[i];
			sum += term;
		}
		return (sum / combinations.size());
	}
	else
	{
		cerr << "This mode not supported!" << endl;
		abort();
	}
}

	//--------------------------------------------------------------------------
	// Perform Bose-Einstein corrections on an event.
	// void shifter::shiftEvent()
	// {
	// 	//auto start = std::chrono::system_clock::now();
	//
	// 	// set up distributions
	// 	int RNG_seed	= (int)paraRdr->getVal("RNG_seed");
	//
	// 	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	// 	default_random_engine generator;
	// 	if ( RNG_seed < 0 )
	// 		generator = default_random_engine (seed);
	// 	else
	// 		generator = default_random_engine (RNG_seed);
	//
	//
	//   // const int n_burn_in = 10000;
	//   const double delta = 0.5;
	//
	//   std::uniform_real_distribution<double> uniform_delta(-delta, delta);
	//   std::uniform_real_distribution<double> uniform(0.0, 1.0);
	//   std::normal_distribution<double>       normal(0.0, 1.0);
	//
	// 	const double R = paraRdr->getVal("RNG_R") / HBARC;
	// 	const double RNG_p0 = paraRdr->getVal("RNG_p0");
	// 	const int number_of_particles = allParticles.size();
	// 	vector<ParticleRecord> allParticles_Original = allParticles;
	//
	// 	constexpr bool check_number_of_shifted_particles = false;
	// 	int iLoop = 0;
	// 	int nLoops = paraRdr->getVal("shifter_nLoops");
	// 	for (iLoop = 0; iLoop < nLoops; iLoop++)
	// 	{
	// 		int number_of_shifted_particles = 0;
	// 		for (int iParticle = 0; iParticle < number_of_particles; iParticle++)
	// 		{
	// 			vector<ParticleRecord> allParticles_copy = allParticles;
	//
	// 			// generate a shifted momentum
	// 			// double x = allParticles.at(iParticle).p.pz();
	// 			double x = allParticles[iParticle].p.pz();
	// 			double y = RNG_p0 * normal(generator);	// corresponds to choice of parameters in random_events.h
	// 			// allParticles_copy.at(iParticle).p.pz(y);
	// 			allParticles_copy[iParticle].p.pz(y);
	//
	// 			// get probability of current configuration
	// 			vector<double> unshifted_pairs = get_pairs( allParticles );
	// 			double P1 = get_probability( R, unshifted_pairs );
	//
	//
	// 			// get probability of shifted configuration
	// 			vector<double> shifted_pairs = get_pairs( allParticles_copy );
	// 			double P2 = get_probability( R, shifted_pairs );
	//
	// 			// choose new configuration (shifted or original)
	// 	    double log_alpha = std::min(0.0, log(P2/P1) + 0.0*(y*y - x*x));
	// 			bool shift_this_particle = (log(uniform(generator)) <= log_alpha);
	// 	    x = shift_this_particle ? y : x;
	//
	// 			if (shift_this_particle) number_of_shifted_particles++;
	//
	// 			// store (possibly new) x to appropriate particle
	// 			// allParticles.at(iParticle).p.pz(x);
	// 			allParticles[iParticle].p.pz(x);
	// 		}
	// 		if ( check_number_of_shifted_particles && number_of_shifted_particles == 0 ) break;
	// 	}
	//
	// 	// cout << "--------------------------------------------------------------------------------" << endl;
	// 	// for (int iParticle = 0; iParticle < number_of_particles; iParticle++)
	// 	// 	cout << "OUT: " << number_of_shifted_events << "   " << iParticle
	// 	// 			<< "   " << allParticles_Original.at(iParticle).p.pz()
	// 	// 			<< "   " << allParticles.at(iParticle).p.pz()
	// 	// 			<< "   " << allParticles_Original.at(iParticle).p.pz()
	// 	// 					- allParticles.at(iParticle).p.pz() << endl;
	// 	vector<double> old_pairs = get_pairs( allParticles_Original );
	// 	vector<double> new_pairs = get_pairs( allParticles );
	//
	// 	for (int iPair = 0; iPair < new_pairs.size(); iPair++)
	// 		cout << "OUT: " << number_of_shifted_events << "   " << iPair
	// 				<< "   " << old_pairs[iPair]
	// 				<< "   " << new_pairs[iPair]
	// 				<< "   " << new_pairs[iPair] - old_pairs[iPair] << "\n";
	//
	// 	cerr << "CHECK: " << standard_deviation( old_pairs )
	// 				<< "   " << standard_deviation( new_pairs ) << "   "
	// 				<< "   " << standard_deviation( new_pairs )
	// 											- standard_deviation( old_pairs ) << "\n";
	//
	// 	cerr << "Finished shifting in " << iLoop << " of " << nLoops << " loops.\n";
	//
	// 	// Done.
	// 	return;
	//
	// }

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
	int RNG_seed	= (int)paraRdr->getVal("RNG_seed");

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator;
	if ( RNG_seed < 0 )
		generator = default_random_engine (seed);
	else
		generator = default_random_engine (RNG_seed);

	std::uniform_real_distribution<double> uniform(0.0, 1.0);
	std::normal_distribution<double>       normal(0.0, 1.0);

	// this RNG is for shuffling particle momenta
	auto rng = std::default_random_engine { std::random_device{}() };

	const double R = paraRdr->getVal("RNG_R") / HBARC;
	const double RNG_p0 = paraRdr->getVal("RNG_p0");
	const int number_of_particles = allParticles.size();
	vector<ParticleRecord> allParticles_Original = allParticles;

	vector<double> current_pairs = get_pairs(allParticles);

	// need comparator for sorting particles by momentum
	auto particleSort = [](const ParticleRecord & p1,
                         const ParticleRecord & p2)
                        { return p1.p.pz() < p2.p.pz(); };

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

			if ( SHIFT_MODE == "CramerVonMises" )
				shift_this_particle = bool( P2 > P1 );

			x = shift_this_particle ? y : x;

// cout << "ratio = " << P2/P1 << ";   " << x << endl;

			if (shift_this_particle)
			{
				current_pairs = shifted_pairs;
				number_of_shifted_particles++;
			}

			// store (possibly new) x to appropriate particle
			allParticles[iParticle].p.pz(x);
		}
		// sort(allParticles.begin(), allParticles.end(), particleSort);
		if ( check_number_of_shifted_particles && number_of_shifted_particles == 0 ) break;
	}

	// std::shuffle(std::begin(allParticles), std::end(allParticles), rng);


// if (true) abort();

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
