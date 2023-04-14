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
#include <unordered_set>

#include "../include/shifter.h"
#include "../include/ParameterReader.h"

using namespace std;

namespace shift_lib
{
	// constexpr const char * SHIFT_MODE = "TRIAL3";

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

		// cout << "SHIFT_MODE = " << SHIFT_MODE << endl;

		// process (i.e., shift and output) this event
		process_event( allParticles_in );

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

		// process (i.e., shift) this event
		process_event( allParticles_in );

		return;
	}

	void shifter::process_event( vector<ParticleRecord> & allParticles_in )
	{
		allParticles = allParticles_in;

		// Perform shifts.
		shiftEvent_efficient();

		// Return shifted results.
		allParticles_in	= allParticles;

		number_of_shifted_events++;

		return;
	}


	void shifter::print(int eventID, vector<ParticleRecord> & allParticles_in, const string & filename)
	{
		ofstream out(filename.c_str(), ios::app);
		out << eventID << "   " << allParticles_in.size() << "\n";
		for (const auto & particle: allParticles_in) out << particle.p.pz() << "\n";
		out.close();
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
			vector<int> combination(K);

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


//------------------------------------------------------------------------------
double shifter::get_probability( const double R, const vector<double> & pair_qzs )
{
	if ( SHIFT_MODE == "Exact" )
	{
		if (pair_qzs.size() != 3)
		{
			cerr << "Cannot run Exact mode with " << pair_qzs.size() << " pairs\n";
			terminate();
		}
		double result = 1.0;
		for (const auto & qz: pair_qzs) result += exp(-0.5*qz*qz*R*R);

		// extra term from Zajc's paper
		result += 2.0*exp(-0.25*R*R*(pair_qzs[0]*pair_qzs[0]
																	+pair_qzs[1]*pair_qzs[1]
																	+pair_qzs[2]*pair_qzs[2]));
		return result;
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
		// use only np-1 independent pairs: first np-1 pairs in list
		const int n = pair_qzs.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		for (int i = 0; i < np-1; i++)
			result *= 1.0 + normalization*exp(-0.5*pair_qzs[i]*pair_qzs[i]*R*R);
		return result;
	}
	//--------------------------------------------------------------------------
	else if ( SHIFT_MODE == "TRIAL2" )
	{
		// use only np-1 independent pairs: 0:1, 1:2, 2:3, ..., np-2:np-1
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
	//--------------------------------------------------------------------------
	else if ( SHIFT_MODE == "TRIAL2b" )
	{
		// use np (np-1 independent) pairs: 0:1, 1:2, 2:3, ..., np-1:np
		const int n = pair_qzs.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		int i = -1;
		for (int i1 = 0; i1 < np - 1; ++i1)
		for (int i2 = i1 + 1; i2 < np; ++i2)
		{
			i++;
			bool include_this_pair = (i2 == i1+1) || (i1 == 0 && i2 == np-1);
			if (!include_this_pair) continue;
			result *= 1.0 + 0.5*(np-1.)*normalization*exp(-0.5*pair_qzs[i]*pair_qzs[i]*R*R);
		}
		return result;
	}
	//--------------------------------------------------------------------------
	else if ( SHIFT_MODE == "TRIAL3" )
	{
		// use only np-1 independent pairs, and cycle over which gets omitted
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
	//--------------------------------------------------------------------------
	else if ( SHIFT_MODE == "TRIAL3b" )
	{
		// use only np-1 independent pairs, and cycle over which gets omitted
		const int n = pair_qzs.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		int i = -1;
		double factor = 0.0;
		vector<double> terms;
		for (int i1 = 0; i1 < np - 1; ++i1)
		for (int i2 = i1 + 1; i2 < np; ++i2)
		{
			i++;
			int di = std::abs(i2-i1);
			bool include_this_pair = (std::min(di, np-di) == 1);
			if (!include_this_pair) continue;
			double term = 1.0 + 0.5*np*normalization*exp(-0.5*pair_qzs[i]*pair_qzs[i]*R*R);
			result *= term;
			factor += 1.0/term;
			terms.push_back(1.0/term);
		}
		for (auto & term: terms) term *= result;
		return *max_element(terms.begin(), terms.end());
	}
	//--------------------------------------------------------------------------
	else if ( SHIFT_MODE == "TRIAL4" )
	{
		// use only np-1 independent pairs, and cycle over all possible combinations
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
	//--------------------------------------------------------------------------
	else if ( SHIFT_MODE == "TRIAL4b" )
	{
		// use only np-1 independent pairs, and cycle over random sample of all possible combinations
		const int n = pair_qzs.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");

		// store all factors for each pair
		vector<double> factors(n);
		for (int i = 0; i < n; i++)
			factors[i] = 1.0 + 0.5*np*normalization*exp(-0.5*pair_qzs[i]*pair_qzs[i]*R*R);

		double sum = 0.0;
		double term_count = 0.0;
		std::string bitmask(np-1, 1); // np-1 leading 1's
		bitmask.resize(n, 0); // n-np+1 trailing 0's

		// fill unordered set with permutations of bitmask
		std::unordered_set<std::string> bitmasks;
		while (bitmasks.size() <= 10*np)
		{
			std::random_shuffle(bitmask.begin(), bitmask.end());
			bitmasks.insert( bitmask );
		}

		// now iterate over permuted bitmasks
		for (const auto & bitmask_permutation: bitmasks)
		{
			double term = 1.0;
			for (int i = 0; i < n; ++i) // [0..n-1] integers
				if (bitmask_permutation[i]) term *= factors[i];
			sum += term;
			term_count++;
		}

// cout << "Exiting here" << endl;

		return (sum / term_count);
	}
	//--------------------------------------------------------------------------
	else if ( SHIFT_MODE == "TRIAL5" )
	{
		// use adjacent particle pairs and next-to-neighbor pairs
		const int n = pair_qzs.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double total = 0.0;
		double normalization = paraRdr->getVal("shifter_norm");
		int maxsep = np/2;
		for (int step = 1; step <= maxsep; step++) // sum over independent pairs (modulo step)
		{
			int i = -1;
			double result = 1.0;
			for (int i1 = 0; i1 < np - 1; ++i1)
			for (int i2 = i1 + 1; i2 < np; ++i2)
			{
				i++;
				int di = std::abs(i2-i1);
				bool include_this_pair = (std::min(di, np-di) == step);
				if (!include_this_pair) continue;
				result *= 1.0 + 0.5*(np-1.0)*normalization*exp(-0.5*pair_qzs[i]*pair_qzs[i]*R*R);
			}
			total += result;
		}
		return total/maxsep;
	}
	//--------------------------------------------------------------------------
	else if ( SHIFT_MODE == "TRIAL5b" )
	{
		// use adjacent particle pairs and next-to-neighbor pairs
		const int n = pair_qzs.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double normalization = paraRdr->getVal("shifter_norm");
		int maxsep = np/2;
		vector<double> totals;
		for (int step = 1; step <= maxsep; step++) // sum over independent pairs (modulo step)
		{
			int i = -1;
			double result = 1.0, factor = 0.0;
			for (int i1 = 0; i1 < np - 1; ++i1)
			for (int i2 = i1 + 1; i2 < np; ++i2)
			{
				i++;
				int di = std::abs(i2-i1);
				bool include_this_pair = (std::min(di, np-di) == step);
				if (!include_this_pair) continue;
				double term = 1.0 + 0.5*np*normalization*exp(-0.5*pair_qzs[i]*pair_qzs[i]*R*R);
				result *= term;
				factor += 1.0/term;
			}
			totals.push_back(factor*result/np);
		}
		// auto rms = [](const vector<double> & v)
    //              { double sum = 0.0;
    //                for (const auto & e: v) sum += e*e;
    //                return sqrt(sum/v.size()); };
		auto mean = [](const vector<double> & v)
                  { double sum = 0.0;
                    for (const auto & e: v) sum += e;
                    return sqrt(sum/v.size()); };
		return mean(totals);
	}
	//--------------------------------------------------------------------------
	else if ( SHIFT_MODE == "RMSscale" )
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
	//--------------------------------------------------------------------------
	else
	{
		cerr << "This mode (" << SHIFT_MODE << ") not supported!" << endl;
		abort();
	}
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
	const double RMSscale = get_RMSscale(allParticles) / HBARC;
// if (true)
// {
// 	cout << "R = " << R << endl;
// 	cout << "RMSscale = " << RMSscale << endl;
// 	terminate();
// }
	const double RNG_p0 = paraRdr->getVal("RNG_p0");
	const int number_of_particles = allParticles.size();
	vector<ParticleRecord> allParticles_Original = allParticles;

	vector<double> current_pairs = get_pairs(allParticles);

	// need comparator for sorting particles by momentum
	auto particleSort = [](const ParticleRecord & p1,
                         const ParticleRecord & p2)
                        { return p1.p.pz() < p2.p.pz(); };

	double proposal_weight = paraRdr->getVal("proposal_weight");
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
			double P1 = (SHIFT_MODE == "RMSscale") ?
									get_probability( RMSscale, current_pairs ) :
									get_probability( R, current_pairs );

			// get probability of shifted configuration
			vector<ParticleRecord> allParticles_with_shift = allParticles;
			allParticles_with_shift[iParticle].p.pz(y);
			vector<double> shifted_pairs = get_shifted_pairs( current_pairs,
																			allParticles_with_shift, iParticle );
			double P2 = (SHIFT_MODE == "RMSscale") ?
									get_probability( RMSscale, shifted_pairs ) :
									get_probability( R, shifted_pairs );

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
		if ( check_number_of_shifted_particles && number_of_shifted_particles == 0 ) break;
	}

	vector<double> old_pairs = get_pairs( allParticles_Original );
	vector<double> new_pairs = get_pairs( allParticles );

	// for (int iPair = 0; iPair < new_pairs.size(); iPair++)
	// 	cout << "OUT: " << number_of_shifted_events << "   " << iPair
	// 			<< "   " << old_pairs[iPair]
	// 			<< "   " << new_pairs[iPair]
	// 			<< "   " << new_pairs[iPair] - old_pairs[iPair] << "\n";
	//
	// cerr << "CHECK: " << standard_deviation( old_pairs )
	// 			<< "   " << standard_deviation( new_pairs ) << "   "
	// 			<< "   " << standard_deviation( new_pairs )
	// 										- standard_deviation( old_pairs ) << "\n";
	//
	// cerr << "Finished shifting in " << iLoop << " of " << nLoops << " loops.\n";

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

	//--------------------------------
	double shifter::get_RMSscale( const vector<ParticleRecord> & particles )
	{
		double result = 0.0;
		const int number_of_particles = particles.size();
		const int number_of_pairs = number_of_particles*(number_of_particles-1)/2;

		for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
		for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
		{
			const double dz = particles[i1].x.pz() - particles[i2].x.pz();
			result += dz*dz;
		}

		return sqrt(0.5*result/number_of_pairs);
	}

}	//End of namespace

//End of file
