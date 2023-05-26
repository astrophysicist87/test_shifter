#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "../include/permanent.h"
#include "../include/shifter.h"
#include "../include/stopwatch.h"
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

	constexpr double TINY       = 1e-6;


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
		for (const auto & particle: allParticles_in)
			out << particle.p.px() << "   " << particle.p.py() << "   " << particle.p.pz() << "\n";
		out.close();
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

	double shifter::permanent(const vector<double> & A, long n) // expects n by n matrix encoded as vector
	{
		// std::cout << "\nA(exact) =\n" << fixed << setprecision(6);
		// for (int i = 0; i < n; ++i)
		// {
		// 	for (int j = 0; j < n; ++j)
		// 		std::cout << "  " << A[i*n+j];
		// 	cout << "\n";
		// }
		// cout << "\n";

    double sum = 0.0;
    double rowsumprod = 0.0, rowsum = 0.0;
    vector<long> chi(n + 1);
    double C = (double)pow((double)2, n);

    // loop all 2^n submatrices of A
    for (long k = 1; k < C; k++)
    {
        rowsumprod = 1.0;
        chi = dec2binarr(k, n); // characteristic vector

        // loop columns of submatrix #k
        for (long m = 0; m < n; m++)
        {
            rowsum = 0.0;

            // loop rows and compute rowsum
            for (long p = 0; p < n; p++)
                rowsum += chi[p] * A[m * n + p];

            // update product of rowsums
            rowsumprod *= rowsum;

            // (optional -- use for sparse matrices)
            if (rowsumprod < TINY) break;
        }

        sum += (double)pow((double)-1, n - chi[n]) * rowsumprod;
    }

    return sum;
}

//------------------------------------------------------------------------------
double shifter::get_probability( const double R, const vector<vector<double>> & qVec,
                                 const vector<double> & BE_distances )
{
	if ( SHIFT_MODE == "HYBRID" )
	{
		double result = get_probability( R, qVec, "FullProduct" );
		double hybrid_cutoff = paraRdr->getVal("hybrid_cutoff");
		if (result > hybrid_cutoff)
			result = get_probability( R, qVec, "Exact" );
		return result;
	}
	else if ( SHIFT_MODE == "Exact" )
	{
		const int n = qVec.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		vector<double> A(np*np);
		{
			int index = 0;
			for (int i = 0; i < np; i++)
			{
				A[i*np+i] = 1.0;
				for (int j = i+1; j < np; j++)
				{
					auto q = qVec[index];
					double q2 = inner_product(q.cbegin(), q.cend(),
																		q.cbegin(),
																		0.0);
					double tmp = exp(-0.25*q2*R*R);
// cout << "Check q: " << q[0] << "  " << q[1] << "  " << q[2] << "  " << q2 << "  " << tmp << "\n";
					if (tmp < TINY) tmp = 0.0;	// make matrix as sparse as possible
					A[i*np+j] = tmp;
					A[j*np+i] = tmp; // matrix is symmetric
					index++;
				}
			}
		}
		return permanent(A, np);
	}
	//--------------------------------------------------------------------------
	// FULL PRODUCT
	else if ( SHIFT_MODE == "FullProduct" )
	{
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		for (const auto & q: qVec)
		{
			double q2 = inner_product(q.cbegin(), q.cend(),
																q.cbegin(),
																0.0);
			result *= 1.0 + normalization*exp(-0.5*q2*R*R);
		}
		return result;
	}
	//--------------------------------------------------------------------------
	// EXPERIMENTAL
	else if ( SHIFT_MODE == "TRIAL3" )
	{
		// use only np-1 independent pairs, and cycle over which gets omitted
		const int n = qVec.size();
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
			auto q = qVec[i];
			double q2 = inner_product(q.cbegin(), q.cend(),
																q.cbegin(),
																0.0);
			double term = 1.0 + 0.5*np*normalization*exp(-0.5*q2*R*R);
			result *= term;
			factor += 1.0/term;
		}
		return factor*result/np;
	}
	//--------------------------------------------------------------------------
	// EXPERIMENTAL
	else if ( SHIFT_MODE == "TRIAL3Fast" )
	{
		// use only np-1 independent pairs, and cycle over which gets omitted
		const int n = qVec.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		double normalization = paraRdr->getVal("shifter_norm");

		auto square = [](double x){return x*x;};
		auto UTindexer = [](int i, int j, int n){return -1 + j - i*(3 + i - 2*n)/2;};
		// auto get_q2 = [](const vector<double> & q)
		// 							{ return inner_product( q.cbegin(), q.cend(), q.cbegin(), 0.0 ); };
		// auto get_q2 = [](const vector<double> & q)
		// 							{ return q[0]*q[0] + q[1]*q[1] + q[2]*q[2]; };

		// initialize with first pair
		// double init = 1.0 + 0.5*np*normalization
		// 										*exp(-0.5*get_q2( qVec[UTindexer(0, np-1, np)] )*R*R);
		double init = 1.0 + 0.5 * np * normalization
                        * square( BE_distances[ UTindexer(0, np-1, np) ] );
		double result = init;
		double factor = 1.0/init;

		// update with other pairs
		for (int i1 = 0; i1 < np - 1; ++i1)
		{
			auto i = UTindexer(i1, i1+1, np);
			double term = 1.0 + 0.5 * np * normalization * square( BE_distances[i] );
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



//------------------------------------------------------------------------------
double shifter::get_probability( const double R, const vector<vector<double>> & qVec, const string & CHOSEN_SHIFT_MODE )
{
	if ( CHOSEN_SHIFT_MODE == "Exact" )
	{
		const int n = qVec.size();
		const int np = static_cast<int>(0.5*(1.0+sqrt(1.0+8.0*n)));
		vector<double> A(np*np);
		{
			int index = 0;
			for (int i = 0; i < np; i++)
			{
				A[i*np+i] = 1.0;
				for (int j = i+1; j < np; j++)
				{
					auto q = qVec[index];
					double q2 = inner_product(q.cbegin(), q.cend(),
																		q.cbegin(),
																		0.0);
					double tmp = exp(-0.25*q2*R*R);
					if (tmp < TINY) tmp = 0.0;	// make matrix as sparse as possible
					A[i*np+j] = tmp;
					A[j*np+i] = tmp; // matrix is symmetric
					index++;
				}
			}
		}
		return permanent(A, np);
	}
	//--------------------------------------------------------------------------
	// FULL PRODUCT
	else if ( CHOSEN_SHIFT_MODE == "FullProduct" )
	{
		double result = 1.0;
		double normalization = paraRdr->getVal("shifter_norm");
		for (const auto & q: qVec)
		{
			double q2 = inner_product(q.cbegin(), q.cend(),
																q.cbegin(),
																0.0);
			result *= 1.0 + normalization*exp(-0.5*q2*R*R);
		}
		return result;
	}
	//--------------------------------------------------------------------------
	else if ( CHOSEN_SHIFT_MODE == "TRIAL3" )
	{
		// use only np-1 independent pairs, and cycle over which gets omitted
		const int n = qVec.size();
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
			auto q = qVec[i];
			double q2 = inner_product(q.cbegin(), q.cend(),
																q.cbegin(),
																0.0);
			double term = 1.0 + 0.5*np*normalization*exp(-0.5*q2*R*R);
			result *= term;
			factor += 1.0/term;
		}
		return factor*result/np;
	}
	//--------------------------------------------------------------------------
	else
	{
		cerr << "This mode (" << CHOSEN_SHIFT_MODE << ") not supported!" << endl;
		abort();
	}
}

//--------------------------------
vector<vector<double>> shifter::get_pairs( const vector<ParticleRecord> & particles )
{
	vector<vector<double>> result;

	const int number_of_particles = particles.size();

	// get all values in vector first
	for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
	for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
	{
		// Vec4 q = particles.at(i1).p - particles.at(i2).p;
		Vec4 q = particles[i1].p - particles[i2].p;
		result.push_back( vector<double>({q.px(), q.py(), q.pz()}) );
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
	// auto rng = std::default_random_engine { std::random_device{}() };

	// avoid replicating previously generated values (cycle through twice to be safe)
	for (const auto & p: allParticles) {normal(generator); normal(generator);}

	const double R = paraRdr->getVal("RNG_R") / HBARC;
	const double RNG_p0 = paraRdr->getVal("RNG_p0");

	const int number_of_particles = allParticles.size();
	vector<ParticleRecord> allParticles_Original = allParticles;

	// vector<vector<double>> current_pairs = get_pairs(allParticles);
	current_pairs = get_pairs(allParticles);
	vector<double> current_BE_distances(current_pairs.size(), 0.0);
	auto get_BE_distance = [R](const vector<double> & q)
	                       { return exp(-0.25 * R * R
                                       * (q[0]*q[0] + q[1]*q[1] + q[2]*q[2])); };
	std::transform( current_pairs.begin(), current_pairs.end(),
                  current_BE_distances.begin(), get_BE_distance );

	shifted_pairs = current_pairs;
	vector<double> shifted_BE_distances = current_BE_distances;

	// get probability of initial configuration
	MatrixPermanent mp(allParticles.size(), 1e-2, true);
	double P1 = (SHIFT_MODE == "AlmostExact") ?
              mp.evaluate(allParticles, current_BE_distances) :
              get_probability( R, current_pairs, current_BE_distances );

	int RNG_xDir 	= paraRdr->getVal("RNG_xDir");
	int RNG_yDir 	= paraRdr->getVal("RNG_yDir");
	int RNG_zDir 	= paraRdr->getVal("RNG_zDir");

	Stopwatch swTotal;

	int iLoop = 0;
	int nLoops = paraRdr->getVal("shifter_nLoops");
	for (iLoop = 0; iLoop < nLoops; iLoop++)	// loop over all particles fixed number of times
	{
		// swTotal.Reset();
		// swTotal.Start();
		cerr << "Loop #" << iLoop << "\n";
		for (int iParticle = 0; iParticle < number_of_particles; iParticle++) // loop over particles, re-sample one at a time
		{
			// generate a shifted momentum
			double x1 = allParticles[iParticle].p.px();
			double y1 = allParticles[iParticle].p.py();
			double z1 = allParticles[iParticle].p.pz();
			double x2 = RNG_xDir ? RNG_p0 * normal(generator) : x1;	// corresponds to choice of parameters in random_events.h
			double y2 = RNG_yDir ? RNG_p0 * normal(generator) : y1;	// corresponds to choice of parameters in random_events.h
			double z2 = RNG_zDir ? RNG_p0 * normal(generator) : z1;	// corresponds to choice of parameters in random_events.h

			// compute shifted configuration
			vector<ParticleRecord> allParticles_with_shift = allParticles;
			allParticles_with_shift[iParticle].p.px(x2);
			allParticles_with_shift[iParticle].p.py(y2);
			allParticles_with_shift[iParticle].p.pz(z2);

			get_shifted_pairs_FAST( shifted_pairs, shifted_BE_distances,
															allParticles_with_shift, iParticle, R );

			// get probability of shifted configuration
			double P2 = (SHIFT_MODE == "AlmostExact") ?
		              mp.evaluate(allParticles_with_shift, shifted_BE_distances) :
		              get_probability( R, shifted_pairs, shifted_BE_distances );

			// for debugging purposes, remove eventually
			if (false && SHIFT_MODE == "Exact")
			{
				#pragma omp critical
				{
					cout << "  Compare unshifted: " << P1
								<< "  " << mp.evaluate(allParticles, current_BE_distances)
								<< "  " << get_probability( R, current_pairs, "Exact" )
								<< "  " << get_probability( R, current_pairs, "FullProduct" ) << "\n";
					cout << "  Compare shifted: " << P2
								<< "  " << mp.evaluate(allParticles_with_shift, shifted_BE_distances)
								<< "  " << get_probability( R, shifted_pairs, "Exact" )
								<< "  " << get_probability( R, shifted_pairs, "FullProduct" ) << "\n";
					// cout << "Seed = " << seed << endl;
				}
				// std::terminate();
			}

			std::cout << "CHECK: " << SHIFT_MODE << "   "
								<< x1 << "   " << y1 << "   " << z1 << "   "
								<< x2 << "   " << y2 << "   " << z2 << "   "
								<< P1 << "   " << P2 << std::endl;
			for (const auto & p: allParticles) std::cout << p.p.pz() << "   ";
			std::cout << std::endl;
			for (const auto & p: allParticles_with_shift) std::cout << p.p.pz() << "   ";
			std::cout << std::endl;
			// std::cout << "Iteration/particle = " << iLoop << " / " << iParticle << ": "
			// 					<< get_probability( R, current_pairs, "Exact" ) << "   "
			// 					<< get_probability( R, shifted_pairs, "Exact" ) << "   "
			// 					<< get_probability( R, current_pairs, "FullProduct" ) << "   "
			// 					<< get_probability( R, shifted_pairs, "FullProduct" ) << "   "
			// 					<< get_probability( R, current_pairs, "TRIAL3" ) << "   "
			// 					<< get_probability( R, shifted_pairs, "TRIAL3" ) << "   " << std::endl;

			// choose new configuration (shifted or original)
			double log_alpha = std::min(0.0, log(P2/P1));
			bool shift_this_particle = (log(uniform(generator)) <= log_alpha);

			// re-set shifted quantities
			if (shift_this_particle)
			{
				// current_pairs = shifted_pairs;		// re-set current configuration
				get_shifted_pairs_FAST( current_pairs, current_BE_distances,
																allParticles_with_shift, iParticle, R );
				P1 = P2;													// re-set probability of current configuration
				allParticles[iParticle].p.px(x2);	// re-set momentum of particle
				allParticles[iParticle].p.py(y2);	// re-set momentum of particle
				allParticles[iParticle].p.pz(z2);	// re-set momentum of particle
			}
			else	// overwrite shifted pairs with current (unshifted)
				get_shifted_pairs_FAST( shifted_pairs, shifted_BE_distances,
																allParticles, iParticle, R );

		}
		// swTotal.Stop();
		// cout << "finished in " << swTotal.printTime() << " s.\n";
	}

	// Done.
	return;
}

		//--------------------------------
		vector<vector<double>> shifter::get_shifted_pairs( const vector<vector<double>> & pairs,
		 																					 const vector<ParticleRecord> & particles,
																							 const int shifted_particle_index )
		{
			vector<vector<double>> result = pairs;
			// auto indexer = [nCols](size_t i, size_t j) { return i*nCols+j; };

			const int number_of_particles = particles.size();

			int iPair = 0;
			for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
			for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
			{
				if (i1==shifted_particle_index || i2==shifted_particle_index)
					result[iPair] = vector<double>({
														particles[i1].p.px() - particles[i2].p.px(),
														particles[i1].p.py() - particles[i2].p.py(),
														particles[i1].p.pz() - particles[i2].p.pz()});
				iPair++;
			}

			return result;
		}

		//--------------------------------
		void shifter::get_shifted_pairs_FAST( vector<vector<double>> & pairs,
																					vector<double> & BE_distances,
																				  const vector<ParticleRecord> & particles,
																			    const int shifted_particle_index,
                                          const double R )
		{
			const int np = particles.size();
			auto & spi = shifted_particle_index;

			// indexes upper-triangular list of pairs given indices (i,j)
			auto UTindexer = [](int i, int j, int n){return -1 + j - i*(3 + i - 2*n)/2;};
			auto get_q2 = [](std::initializer_list<double> q)
										{ return inner_product( q.begin(), q.end(), q.begin(), 0.0 ); };

			for (int i1 = 0; i1 < spi; ++i1)
			{
				int iPair = UTindexer(i1, spi, np);
				auto q = { particles[i1].p.px() - particles[spi].p.px(),
									 particles[i1].p.py() - particles[spi].p.py(),
									 particles[i1].p.pz() - particles[spi].p.pz() };
				BE_distances[iPair] = exp(-0.25 * R * R * get_q2(q));
				pairs[iPair].assign(q);
			}

			for (int i1 = spi+1; i1 < np; ++i1)
			{
				int iPair = UTindexer(spi, i1, np);
				auto q = { particles[spi].p.px() - particles[i1].p.px(),
									 particles[spi].p.py() - particles[i1].p.py(),
									 particles[spi].p.pz() - particles[i1].p.pz() };
				BE_distances[iPair] = exp(-0.25 * R * R * get_q2(q));
				pairs[iPair].assign(q);
			}
			return;
		}

}	//End of namespace

//End of file
