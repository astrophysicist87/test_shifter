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

constexpr bool BE_VERBOSE = false;

constexpr double MM2FM = 1.e12;
constexpr double FM2MM = 1.e-12;
constexpr double HBARC = 0.19733;

constexpr double COMPRELERR = 1e-10;
constexpr double COMPFACMAX = 1000.;
constexpr int    NCOMPSTEP  = 10;


void shifter::initialize_all( ParameterReader * paraRdr_in,
	const vector<ParticleRecord> & allParticles_in )
{
	include_pair_density = true;

	// Load parameters
	paraRdr			= paraRdr_in;

	// Load particle
	allParticles	= allParticles_in;

	// Perform shifts
	shiftEvent();

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



//--------------------------------------------------------------------------
// Perform Bose-Einstein corrections on an event.

void shifter::shiftEvent()
{
	// Reset list of identical particles.
	//allParticles.resize(0);

	//auto start = std::chrono::system_clock::now();

	// Loop through event record to store copies of current species.
	/*for (int i = 0; i < event.size(); ++i)
		if ( event[i].id() == idNow and (event[i].isFinal() or debugging) )
		{
			allParticles.push_back(
			shifterHadron( idNow, i, event[i].p(), event[i].m(), event[i].vProd() ) );
		}
	*/

	// Loop through pairs of identical particles and find shifts.
	bool enoughPairsToProceed = setSortedPairs( allParticles );
	if ( enoughPairsToProceed )
	{
		cout << "shifterCheck: NPair = " << sortedPairs.size() << endl;
		shiftPairs_mode1();
	}

	int iPair = 0;
	for (const auto & thisPair: pairs_sorted_by_abs_qz)
		cout << "Unshifted: " << iPair++ << "   " << thisPair.first << "   "
				<< thisPair.second.first << "   " << thisPair.second.second << endl;



	// Add in shifts without compensations
	allParticles_Shifted = allParticles;
	for (auto & thisParticle: allParticles_Shifted)
		thisParticle.p += thisParticle.pShift;

	// Reconstruct original vs. new qz distributions
	enoughPairsToProceed = setSortedPairs( allParticles_Shifted );
	iPair = 0;
	for (const auto & thisPair: pairs_sorted_by_abs_qz)
		cout << "Shifted: " << iPair++ << "   " << thisPair.first << "   "
				<< thisPair.second.first << "   " << thisPair.second.second << endl;


	// Must have at least two pairs to carry out compensation.
	const int nParticles = allParticles.size();
	if (nParticles < 2) return;

	// Add in compensations until energy is conserved
	double eSumOriginal = 0.;
	double eSumShifted  = 0.;
	double eDiffByComp  = 0.;
	for (auto & particle: allParticles)
	{
		eSumOriginal  += particle.p.e();
		particle.p    += particle.pShift;
		particle.p.e( sqrt( particle.p.pAbs2() + particle.m2 ) );
		eSumShifted   += particle.p.e();
		eDiffByComp   += dot3( particle.pComp, particle.p) / particle.p.e();
	}
	


	constexpr bool perform_compensation = true;

	// Iterate compensation shift until convergence.
	int iStep = 0;
	while ( perform_compensation
	&& abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal
	&& abs(eSumShifted - eSumOriginal) < COMPFACMAX * abs(eDiffByComp)
	&& iStep < NCOMPSTEP )
	{
		++iStep;
		double compFac   = (eSumOriginal - eSumShifted) / eDiffByComp;
		eSumShifted      = 0.;
		eDiffByComp      = 0.;
		for (auto & particle: allParticles)
		{
			particle.p    += compFac * particle.pComp;
			particle.p.e( sqrt( particle.p.pAbs2() + particle.m2 ) );
			eSumShifted   += particle.p.e();
			eDiffByComp   += dot3( particle.pComp, particle.p) / particle.p.e();
		}
	}


	// Output final qz distribution
	enoughPairsToProceed = setSortedPairs( allParticles );
	iPair = 0;
	for (const auto & thisPair: pairs_sorted_by_abs_qz)
		cout << "Compensated: " << iPair++ << "   " << thisPair.first << "   "
				<< thisPair.second.first << "   " << thisPair.second.second << endl;



	constexpr bool check_for_bad_events = false;

	// Error if no convergence, and then return without doing BE shift.
	// However, not grave enough to kill event, so return true.
	if ( perform_compensation
	and check_for_bad_events
	and abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal )
	{
		//infoPtr->errorMsg("Warning in shifter::shiftEvent: no consistent BE shift topology found, so skip BE");
		cout << setprecision(16) << "shifterCheck: This event did not pass! Check: "
				<< abs(eSumShifted - eSumOriginal) << " < " << COMPRELERR * eSumOriginal << "\n";
		return;
	}
	else
	{
		cout << setprecision(16) << "shifterCheck: This event passes! Check: "
				<< abs(eSumShifted - eSumOriginal) << " < " << COMPRELERR * eSumOriginal << "\n";
	}


	// Store new particle copies with shifted momenta.
	/*for (int i = 0; i < nStored[9]; ++i)
	{
		int iNew = event.copy( allParticles.at(i).iPos, 99);
		event[ iNew ].p( allParticles.at(i).p );
	}*/

	//auto end = std::chrono::system_clock::now();

	//std::chrono::duration<double> elapsed_seconds = end-start;
	//std::cout << "shifterCheck: elapsed time: " << elapsed_seconds.count() << " s" << "\n";

	// Done.
	return;

}

//--------------------------------
bool shifter::setSortedPairs( const vector<ParticleRecord> & particles_to_sort )
{
	// Reset.
	sortedPairs.clear();

	pairs_sorted_by_qzPRF.clear();
	pairs_sorted_by_abs_qzPRF.clear();

	pairs_sorted_by_qz.clear();
	pairs_sorted_by_abs_qz.clear();

	const int number_of_particles = particles_to_sort.size();

	// get all values in vector first
	for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
	for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
	{
		Vec4 q = particles_to_sort.at(i1).p - particles_to_sort.at(i2).p;
		Vec4 qPRF = q;
		qPRF.bstback( 0.5*( particles_to_sort.at(i1).p + particles_to_sort.at(i2).p ) );

		/*sortedPairs.push_back(
			std::make_pair(
				sqrt( m2(particles_to_sort.at(i1).p, particles_to_sort.at(i2).p) - m2Pair[iTab] ),
				std::make_pair(i1, i2)
			)
		);*/
		sortedPairs.push_back( std::make_pair( qPRF.pAbs(), std::make_pair(i1, i2) ) );
		pairs_sorted_by_qzPRF.push_back( std::make_pair( qPRF.pz(), std::make_pair(i1, i2) ) );
		pairs_sorted_by_abs_qzPRF.push_back( std::make_pair( abs(qPRF.pz()), std::make_pair(i1, i2) ) );

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
	pairs_sorted_by_abs_qzPRF.push_back( std::make_pair( 0.0, std::make_pair(-1, -1) ) );
	pairs_sorted_by_abs_qz.push_back( std::make_pair( 0.0, std::make_pair(-1, -1) ) );

	// THEN sort them (sorts on first column in ascending order automatically)
	sort( sortedPairs.begin(), sortedPairs.end() );
	sort( pairs_sorted_by_qzPRF.begin(), pairs_sorted_by_qzPRF.end() );
	sort( pairs_sorted_by_abs_qzPRF.begin(), pairs_sorted_by_abs_qzPRF.end() );
	sort( pairs_sorted_by_qz.begin(), pairs_sorted_by_qz.end() );
	sort( pairs_sorted_by_abs_qz.begin(), pairs_sorted_by_abs_qz.end() );

	return (true);
}


void shifter::shiftPairs_mode1()
{
	vector< pair< double, double > > LHS, RHS, RHS_derivatives;
	evaluate_shift_relation_at_pair( pairs_sorted_by_abs_qz, LHS, RHS, RHS_derivatives );


	const double dqz = 0.01;	// GeV
	for (double qzVal = 0.0; qzVal < 20.0 + 0.5*dqz; qzVal+=dqz)
		cout << "Effective source: " << dqz << "   "
				<< evaluate_effective_source( pairs_sorted_by_abs_qz, qzVal ) << endl;

if (1) exit(8);

	compute_shifts( pairs_sorted_by_abs_qz, LHS, RHS, RHS_derivatives );

	/*
	// some output to check stuff
	const int npairs = LHS.size();
	cout << "sizes: "
			<< LHS.size() << "   "
			<< RHS.size() << "   "
			<< pairs_sorted_by_abs_qz.size() << "   "
			//<< pairShifts.size()
			<< endl;
	for (int i = 1; i < npairs; i++)
	{
		const auto & thisPair = pairs_sorted_by_abs_qz.at(i);
		const double this_qz = thisPair.first;
		const int this1 = thisPair.second.first;
		const int this2 = thisPair.second.second;
		const double thisPair_shift = pairShifts.at(i-1);

		cout << setprecision(16) << "CHECK: "
				<< LHS.at(i).first << "   "
				<< LHS.at(i).second << "   "
				<< RHS.at(i).second << "   "
				<< thisPair_shift
				<< endl;
	}

	if (1) exit (8);
	*/

	return;
}





//-------------------------------------
// Compute the unshifted pair integrals
// at each pair's relative momentum.
void shifter::evaluate_shift_relation_at_pair(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			vector< pair< double, double > > & LHS,
			vector< pair< double, double > > & RHS,
			vector< pair< double, double > > & RHS_derivatives
			)
{

	//-------
	// Reset.
	LHS.clear();
	RHS.clear();
	RHS_derivatives.clear();

	// Make sure pair density is stored, if needed.
	set_pair_density( sorted_list_of_pairs );

	//------------------
	// Set LHS integral.
	set_LHS( sorted_list_of_pairs, LHS );

	//------------------
	// Set RHS integral.
	set_RHS( sorted_list_of_pairs, RHS, RHS_derivatives );

	return;
}



void shifter::set_pair_density(	const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs )
{
	// Reset.
	denBar.resize(0);

	if ( include_pair_density )
	{
		// for counting number of particles below given |q_z|
		for (int iPair = 0; iPair < (int)sorted_list_of_pairs.size()-1; ++iPair)
			denBar.push_back( 0.5 / ( sorted_list_of_pairs.at(iPair+1).first - sorted_list_of_pairs.at(iPair).first ) );
	}
	else
	{
		for (int iPair = 0; iPair < (int)sorted_list_of_pairs.size()-1; ++iPair)
			denBar.push_back( 1.0 );
	}
	

	return;
}


//-------------------------------------
// Set lefthand side of shift relation.
void shifter::set_LHS(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			vector< pair< double, double > > & LHS )
{

	// LHS is symmetric integral about q = 0 up to q = q_pair, for each pair
	// ==>> most efficient to start with smallest |q_pair| and work my way up
	const int npairs = sorted_list_of_pairs.size();
	LHS.reserve( npairs );

	for (const auto & thisPair : sorted_list_of_pairs)
	{
		const double LHS_integral = evaluate_LHS( sorted_list_of_pairs, thisPair.first );
		LHS.push_back( std::make_pair( thisPair.first, LHS_integral ) );
		cout << "CHECK LHS: " << thisPair.first << "   " << LHS_integral << endl;
	}

	return;
}

double shifter::evaluate_LHS(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			double qz_in )
{
	int    pairIndex    = 1;
	double previous_qz  = 0.0;
	double this_qz      = 0.0;
	double LHS_integral = 0.0;

	double qz = abs(qz_in);
	if ( qz < 1.e-20 ) return (0.0);

	while ( this_qz < qz )
	{
		//cout << "checkpoint: " << previous_qz << "   " << this_qz << endl;
		previous_qz   = sorted_list_of_pairs.at(pairIndex-1).first;
		this_qz       = sorted_list_of_pairs.at(pairIndex).first;
		LHS_integral += 2.0 * ( this_qz - previous_qz ) * denBar.at(pairIndex-1);
		pairIndex++;
	}

	//LHS_integral += 2.0 * ( qz - previous_qz ) * denBar.at(pairIndex-1);

	return (LHS_integral);
}


//-------------------------------------
// Set righthand side of shift relation.
void shifter::set_RHS(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			vector< pair< double, double > > & RHS,
			vector< pair< double, double > > & RHS_derivatives )
{

	// RHS
	const int npairs = sorted_list_of_pairs.size();
	RHS.reserve( npairs );

	for (const auto & thisPair : sorted_list_of_pairs)
	{
		double RHS_derivative = 0.0;
		const double RHS_integral = evaluate_RHS( sorted_list_of_pairs, RHS, thisPair, thisPair.first, RHS_derivative );
		RHS.push_back( std::make_pair( thisPair.first, RHS_integral ) );
		RHS_derivatives.push_back( std::make_pair( thisPair.first, RHS_derivative ) );
		cout << "CHECK RHS: " << thisPair.first << "   " << RHS_integral << endl;
	}


/*cout << "Symmetry check of RHS: " << endl;
const double delQ = 2.0*sorted_list_of_pairs.back().first / 1000.0;
double tmp = 0.0, tmp_d = 0.0;
for (int i = -500; i <= 500; i++)
{
	tmp = delQ * i;
	cout << "fineRHS: " << tmp << "   "
			<< evaluate_RHS( sorted_list_of_pairs, RHS,
								sorted_list_of_pairs.front(),
								tmp, tmp_d ) << endl;
}
if (1) exit(8);*/


	return;
}

double shifter::evaluate_RHS(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			const vector< pair< double, double > > & RHS,
			const pair< double, pair <int,int> > & thisPair,
			const double qz_in, double & RHS_derivative )
{
	int    pairIndex = 1;
	double lower_qz  = sorted_list_of_pairs.at(pairIndex-1).first;
	double upper_qz  = sorted_list_of_pairs.at(pairIndex).first;

	//const double qz = thisPair.first;

	double qz = abs(qz_in);
	if ( qz < 1.e-20 ) return (0.0);

	const int this1 = thisPair.second.first;
	const int this2 = thisPair.second.second;

	// recycle previously computed RHS to save time
	while ( upper_qz < qz and pairIndex < sorted_list_of_pairs.size() - 1 )
	{
		lower_qz = sorted_list_of_pairs.at(pairIndex).first;
		upper_qz = sorted_list_of_pairs.at(++pairIndex).first;
	}

	lower_qz = sorted_list_of_pairs.at(pairIndex-1).first;
	upper_qz = sorted_list_of_pairs.at(pairIndex).first;
	double RHS_integral = RHS.at(pairIndex-1).second;

	// the constant piece (integrand assumed to be symmetric)
	//RHS_integral += 2.0 * ( upper_qz - lower_qz ) * denBar.at(pairIndex-1);
	RHS_integral += 2.0 * ( qz - lower_qz ) * denBar.at(pairIndex-1);
	RHS_derivative = 2.0 * denBar.at(pairIndex-1);

	int npairs_in_average = 0;
	double RHS_BE_enhancement = 0.0;
	double RHS_BE_enhancement_derivative = 0.0;

	// the BE enhancement piece
	for (const auto & iPair : sorted_list_of_pairs)
	{
		const int i1 = iPair.second.first;
		const int i2 = iPair.second.second;

		if ( i1<0 or i2<0 ) continue;
		//if ( this1 != i1 and this2 != i2 ) continue;
		//if ( this1 != i1 or  this2 != i2 ) continue;

		Vec4 xDiff = ( allParticles.at(i1).x - allParticles.at(i2).x ) / HBARC;

		const double Delta_z = xDiff.pz();

		RHS_BE_enhancement += 2.0 * ( sin(qz*Delta_z) - sin(lower_qz*Delta_z) )
							* denBar.at(pairIndex-1) / Delta_z;

		RHS_BE_enhancement_derivative += 2.0 * cos(qz*Delta_z) * denBar.at(pairIndex-1);

		npairs_in_average++;
	}

	RHS_BE_enhancement /= npairs_in_average;
	RHS_BE_enhancement_derivative /= npairs_in_average;
	RHS_integral += RHS_BE_enhancement;
	RHS_derivative += RHS_BE_enhancement_derivative;
	pairIndex++;

	return (RHS_integral);
}



/*double shifter::Newtons_Method( const double a, const double b )
{
	const double ACCURACY = 1.e-6;
	const int MAXTRIES = 100;

	const double guess1 = -1.0/b;
	const double guess2 = 1.0/b;

	const double f1 = guess1*b + sin(b*(a+guess1));
	const double f2 = guess2*b + sin(b*(a+guess2));

	//cout << f1 << "   " << f2 << endl;

	// Solve equation given by x*b + sin((a+x)*b) == 0
	const double initial_guess = ( f1*f1 < f2*f2 ) ? guess1 : guess2;

	double x = initial_guess;
	double f = x*b + sin(b*(a+x));

	int ntries = 0;
	while ( abs(f) > ACCURACY and ntries < MAXTRIES )
	{
		double fp = b*(1.0 + cos(b*(a+x)));
		//cout << setprecision(24) << "ntries = " << ntries << ": " << x << "   " << f << "   " << fp << endl;
		if ( abs(fp) < 1.e-100 ) break;
		f = x*b + sin(b*(a+x));
		x -= f / fp;
		ntries++;
	}


	if ( ntries == MAXTRIES ) cout << "WARNING: maximum number of tries reached! a=" << a << ", b=" << b << "; root x = " << x << endl;


	return (x);
}*/


double shifter::compute_shift(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			const vector< pair< double, double > > & LHS,
			const vector< pair< double, double > > & RHS,
			const vector< pair< double, double > > & RHS_derivatives,
			int iPair )
{
	const double ACCURACY = 1.e-6;
	const int MAXTRIES = 10;

	const double qz0 = LHS.at(iPair).first;
	const auto & thisPair = sorted_list_of_pairs.at(iPair);
	const double LHS_thisPair = LHS.at(iPair).second;

	// Solve equation given by LHS(qz0) - RHS(qz0 + x) == 0
	const double initial_guess = 0.0;

	/*cout << setprecision(16) << "Check shift computation: "
			<< iPair << "   " << qz0 << "   "
			<< LHS.at(iPair).second << "   "
			<< RHS.at(iPair).second << endl;*/

	double x = initial_guess;
	double f = RHS.at(iPair).second - LHS_thisPair;

	double fp = RHS_derivatives.at(iPair).second;
	int ntries = 0;
	//cout << setprecision(16) << "ntries = " << ntries << ": " << x << "   " << f << "   " << fp << endl;
	while ( abs(f) > ACCURACY and ntries < MAXTRIES )
	{
		if ( abs(fp) < 1.e-100 ) break;
		x -= f / fp;
		f = evaluate_RHS( sorted_list_of_pairs, RHS, thisPair, qz0 + x, fp ) - LHS_thisPair;
		ntries++;
		//cout << setprecision(16) << "ntries = " << ntries << ": " << x << "   " << f << "   " << fp << endl;
	}

	// if we haven't converged yet, try bisection instead
	if ( ntries == MAXTRIES )
	{
		//cout << "WARNING: maximum number of tries reached! Q=" << qz0 << ": LHS="
		//		<< LHS.at(iPair).second << ", RHS=" << RHS.at(iPair).second << "; root x = " << x << endl;

		const int NMAX = 10*MAXTRIES;
		const double TOL = ACCURACY;
		int N=1;
		double fpdummy = 0.0;
		double a = LHS.front().first - qz0, b = LHS.back().first - qz0;
		double sa = sgn( evaluate_RHS( sorted_list_of_pairs, RHS, thisPair, 0.0, fpdummy ) - LHS_thisPair );
		double sb = sgn( evaluate_RHS( sorted_list_of_pairs, RHS, thisPair, b, fpdummy ) - LHS_thisPair );
		double c = 0.0;
		double fnew = 0.0;
		while (N <= NMAX) // limit iterations to prevent infinite loop
		{
			//cout << "--> root in between " << a << " and " << b << endl;

			fnew = evaluate_RHS( sorted_list_of_pairs, RHS, thisPair, qz0 + c, fpdummy ) - LHS_thisPair;
			if ( abs(fnew) < TOL and 0.5*(b - a) < TOL )
			{
				x = c;
				break;
			}
			N++;
			if ( sgn(fnew) == sa )
			{
				a = c;
				sa = sgn( evaluate_RHS( sorted_list_of_pairs, RHS, thisPair, qz0 + a, fpdummy ) - LHS_thisPair );
			}
			else
			{
				b = c; // new interval
				sb = sgn( evaluate_RHS( sorted_list_of_pairs, RHS, thisPair, qz0 + b, fpdummy ) - LHS_thisPair );
			}
			c=0.5*(a + b); // new midpoint
		}
		
	}






//cout << "Made it to line = " << __LINE__ << endl;
/*const int i1 = thisPair.second.first;
const int i2 = thisPair.second.second;
Vec4 xDiff = ( allParticles.at(i1).x - allParticles.at(i2).x ) / HBARC;
const double Delta_z = xDiff.pz();
const double check_root = Newtons_Method(qz0, Delta_z);
cout << "Newton's method gives " << check_root << endl;*/

cout << thisPair.first << "   " << x << endl;

	return (x);

}



void shifter::compute_shifts(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			const vector< pair< double, double > > & LHS,
			const vector< pair< double, double > > & RHS,
			const vector< pair< double, double > > & RHS_derivatives
			)
{
	const int npairs = sorted_list_of_pairs.size();

	pairShifts.clear();
	this_pair_shifted.clear();

	pairShifts.reserve( npairs );
	this_pair_shifted.reserve( npairs );

	// Compute the pair shifts themselves (skip i=0 case, not physical pair).
	for (int i = 1; i < npairs; i++)
	{
		const auto & thisPair = sorted_list_of_pairs.at(i);
		const double this_qz = thisPair.first;
		const int this1 = thisPair.second.first;
		const int this2 = thisPair.second.second;

		Vec4 xDiff = ( allParticles.at(this1).x - allParticles.at(this2).x ) / HBARC;
		const double Delta_z = xDiff.pz();	// units are 1/GeV here

		const double thisPair_shift
						= compute_shift( sorted_list_of_pairs, LHS, RHS, RHS_derivatives, i );

		pairShifts.push_back( thisPair_shift );
		this_pair_shifted.push_back( true );
	}

	// Store them for each pair.
	int pairIndex = 0;
	int number_of_pairs_shifted = 0, number_of_pairs_not_shifted = 0;
	for (const auto & iPair : sorted_list_of_pairs)
	{
		//const double this_qz = iPair.first;
		const int i1 = iPair.second.first;
		const int i2 = iPair.second.second;

		// skip first unphysical "pair"
		if (i1<0 or i2<0) continue;
		const double this_qz = allParticles.at(i1).p.pz() - allParticles.at(i2).p.pz();


//cout << "CHECK SIGNS: " << this_qz << "   " << allParticles.at(i1).p.pz() - allParticles.at(i2).p.pz() << endl;
		constexpr bool rescale_pair_momenta = true;

		const double net_qz_shift = pairShifts.at(pairIndex);
		const double factor = 0.5 * net_qz_shift / this_qz;

		Vec4 pDiff = factor * (allParticles.at(i1).p - allParticles.at(i2).p);

		if ( rescale_pair_momenta or this_pair_shifted.at(pairIndex) )
		{
			// Compute appropriate shift for pair
			number_of_pairs_shifted++;

			allParticles.at(i1).pShift += pDiff;
			allParticles.at(i2).pShift -= pDiff;

			if ( rescale_pair_momenta )
			{
				// add symmetrically to both momenta
				pDiff = 0.5*(allParticles.at(i1).p + allParticles.at(i2).p);
				allParticles.at(i1).pComp += pDiff;
				allParticles.at(i2).pComp += pDiff;
			}
		}
		else
		{
			// Use computed shift for compensation
			number_of_pairs_not_shifted++;

			//*/
			//Vec4 pDiff = allParticles.at(i1).p - allParticles.at(i2).p;
			allParticles.at(i1).pComp += pDiff;
			allParticles.at(i2).pComp -= pDiff;
		}

		pairIndex++;

	}

	if ( BE_VERBOSE or true )
	{
		cout << "number_of_pairs_shifted = " << number_of_pairs_shifted << endl;
		cout << "number_of_pairs_not_shifted = " << number_of_pairs_not_shifted << endl;
	}

	return;
}


// current effective source is 1 + < cos(qz Delta_z) >
double shifter::evaluate_effective_source(
				const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
				const double qz )
{
	const int npairs_in_average = 0;
	double effective_source = 0.0;

	// the BE enhancement piece
	for (const auto & iPair : sorted_list_of_pairs)
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



//End of file
