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
	bool enoughPairsToProceed = setSortedPairs();
	if ( enoughPairsToProceed )
	{
		cout << "shifterCheck: NPair = " << sortedPairs.size() << endl;
		shiftPairs_mode1();
	}

	/*
	// Must have at least two pairs to carry out compensation.
	//if (nStored[9] < 2) return true;

	// Shift momenta and recalculate energies.
	double eSumOriginal = 0.;
	double eSumShifted  = 0.;
	double eDiffByComp  = 0.;
	for (int i = 0; i < nStored[9]; ++i)
	{
		eSumOriginal  += allParticles.at(i).p.e();
		allParticles.at(i).p += allParticles.at(i).pShift;
		allParticles.at(i).p.e( sqrt( allParticles.at(i).p.pAbs2() + allParticles.at(i).m2 ) );
		eSumShifted   += allParticles.at(i).p.e();
		eDiffByComp   += dot3( allParticles.at(i).pComp, allParticles.at(i).p) / allParticles.at(i).p.e();
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
		for (int i = 0; i < nStored[9]; ++i)
		{
			allParticles.at(i).p += compFac * allParticles.at(i).pComp;
			allParticles.at(i).p.e( sqrt( allParticles.at(i).p.pAbs2() + allParticles.at(i).m2 ) );
			eSumShifted   += allParticles.at(i).p.e();
			eDiffByComp   += dot3( allParticles.at(i).pComp, allParticles.at(i).p) / allParticles.at(i).p.e();
		}
	}

	constexpr bool check_for_bad_events = false;

	// Error if no convergence, and then return without doing BE shift.
	// However, not grave enough to kill event, so return true.
	if ( perform_compensation
	and check_for_bad_events
	and abs(eSumShifted - eSumOriginal) > COMPRELERR * eSumOriginal )
	{
		infoPtr->errorMsg("Warning in shifter::shiftEvent: no consistent BE shift topology found, so skip BE");
		cout << setprecision(16) << "shifterCheck: This event did not pass! Check: "
			<< abs(eSumShifted - eSumOriginal) << " < " << COMPRELERR * eSumOriginal << "\n";
		return true;
	}
	else
	{
		cout << setprecision(16) << "shifterCheck: This event passes! Check: "
		<< abs(eSumShifted - eSumOriginal) << " < " << COMPRELERR * eSumOriginal << "\n";
	}


	// Store new particle copies with shifted momenta.
	for (int i = 0; i < nStored[9]; ++i)
	{
		int iNew = event.copy( allParticles.at(i).iPos, 99);
		event[ iNew ].p( allParticles.at(i).p );
	}
	*/

	//auto end = std::chrono::system_clock::now();

	//std::chrono::duration<double> elapsed_seconds = end-start;
	//std::cout << "shifterCheck: elapsed time: " << elapsed_seconds.count() << " s" << "\n";

	// Done.
	return;

}

//--------------------------------
bool shifter::setSortedPairs( )
{
	// Reset.
	sortedPairs.clear();

	pairs_sorted_by_qzPRF.clear();
	pairs_sorted_by_abs_qzPRF.clear();

	pairs_sorted_by_qz.clear();
	pairs_sorted_by_abs_qz.clear();

	const int number_of_particles = allParticles.size();

	// get all values in vector first
	for (int i1 = 0; i1 < number_of_particles - 1; ++i1)
	for (int i2 = i1 + 1; i2 < number_of_particles; ++i2)
	{
		Vec4 q = allParticles.at(i1).p - allParticles.at(i2).p;
		Vec4 qPRF = q;
		qPRF.bstback( 0.5*( allParticles.at(i1).p + allParticles.at(i2).p ) );

		/*sortedPairs.push_back(
			std::make_pair(
				sqrt( m2(allParticles.at(i1).p, allParticles.at(i2).p) - m2Pair[iTab] ),
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
	vector< pair< double, double > > LHS, RHS;
	evaluate_shift_relation_at_pair( pairs_sorted_by_abs_qz, LHS, RHS );

	// some output to check stuff
	const int npairs = LHS.size();
	cout << "sizes: " << LHS.size() << "   " << pairs_sorted_by_abs_qz.size() << endl;
	for (int i = 1; i < npairs; i++)
	{
		const auto & thisPair = pairs_sorted_by_abs_qz.at(i);
		const double this_qz = thisPair.first;
		const int this1 = thisPair.second.first;
		const int this2 = thisPair.second.second;
		//Vec4 xDiff = ( allParticles.at(this1).x - allParticles.at(this2).x ) / HBARC;
		//const double Delta_z = xDiff.pz();	// units are 1/GeV here

		//const double thisPair_shift = Newtons_Method( this_qz, Delta_z );

		cout << setprecision(24) << "CHECK: "
				<< LHS.at(i).first << "   "
				<< LHS.at(i).second << "   "
				<< RHS.at(i).second << "   "
				//<< thisPair_shift
				<< endl;
	}

	if (1) exit(8);

	compute_shifts( pairs_sorted_by_abs_qz );

	return;
}



void shifter::compute_shifts(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs
			)
{
	const int npairs = sorted_list_of_pairs.size();

	vector<double> pairShifts;
	vector<bool> this_pair_shifted;

	pairShifts.reserve( npairs );
	this_pair_shifted.reserve( npairs );

	// Compute the pair shifts themselves.
	for (int i = 0; i < npairs; i++)
	{
		const auto & thisPair = sorted_list_of_pairs.at(i);
		const double this_qz = thisPair.first;
		const int this1 = thisPair.second.first;
		const int this2 = thisPair.second.second;
		Vec4 xDiff = ( allParticles.at(this1).x - allParticles.at(this2).x ) / HBARC;
		const double Delta_z = xDiff.pz();	// units are 1/GeV here

		const double thisPair_shift = Newtons_Method( this_qz, Delta_z );
		pairShifts.push_back( thisPair_shift );
		if (this1==0 or this2==0)	// choose a particle to track
			cout << setprecision(24) << "CHECK: "
					<< this1 << "   " << this2 << "   " << this_qz << "   "
					<< Delta_z*HBARC << "   " << thisPair_shift << endl;
		this_pair_shifted.push_back( true );
	}

	// Store them for each pair.
	int pairIndex = 0;
	int number_of_pairs_shifted = 0, number_of_pairs_not_shifted = 0;
	for (const auto & iPair : sorted_list_of_pairs)
	{
		const double this_qz = iPair.first;
		const int i1 = iPair.second.first;
		const int i2 = iPair.second.second;

		constexpr bool rescale_pair_momenta = true;

		const double net_qz_shift = 0.5*pairShifts.at(pairIndex);
		const double factor = net_qz_shift / this_qz;

		// Add shifts to sum. (Energy component dummy.)
		//Vec4 pDiff(0.0, 0.0, 0.0, net_qz_shift);
		Vec4 pDiff = factor * (allParticles.at(i1).p - allParticles.at(i2).p);
		if (i1==0 or i2==0)	// choose a particle to track
			cout << setprecision(24) << "CHECK: "
					<< i1 << "   " << i2 << "   " << this_qz << "   "
					<< net_qz_shift << "   " << factor << "   " << pDiff;

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




//-------------------------------------
// Compute the unshifted pair integrals
// at each pair's relative momentum.
void shifter::evaluate_shift_relation_at_pair(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			vector< pair< double, double > > & LHS,
			vector< pair< double, double > > & RHS
			)
{

	//-------
	// Reset.
	LHS.clear();
	RHS.clear();

	// Make sure pair density is stored, if needed.
	set_pair_density( sorted_list_of_pairs );

	//------------------
	// Set LHS integral.
	set_LHS( sorted_list_of_pairs, LHS );

	//------------------
	// Set RHS integral.
	set_RHS( sorted_list_of_pairs, RHS );

	return;
}

double shifter::Newtons_Method( const double a, const double b )
{
	const double ACCURACY = 1.e-6;
	const int MAXTRIES = 100;

	/*const double guess1 = -1.0/b;
	const double guess2 = 0.0;
	const double guess3 = 1.0/b;

	const double f1 = guess1*b + sin(b*(a+guess1));
	const double f2 = guess2*b + sin(b*(a+guess2));
	const double f3 = guess3*b + sin(b*(a+guess3));

	bool f1_lt_f2 = (f1*f1 < f2*f2);
	bool f1_lt_f3 = (f1*f1 < f3*f3);
	bool f2_lt_f3 = (f2*f2 < f3*f3);

	// Solve equation given by x*b + sin((a+x)*b) == 0
	const double initial_guess = guess1 * static_cast<double>(  f1_lt_f2 and f1_lt_f3 )
									+ guess2 * static_cast<double>( not f1_lt_f2 and f2_lt_f3 )
									+ guess3 * static_cast<double>( not f1_lt_f3 and not f2_lt_f3 );*/

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

	//if ( abs(x) > 0.1*a ) cout << setprecision(24) << "WARNING: large shift! a=" << a << ", b=" << b << "; root x = " << x << endl;

	return (x);
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
	
	// Density must go to zero eventually.
	//denBar.back() = 0.0;
	//denBar.push_back( 0.0 );

	/*cout << "Sizes: " << sorted_list_of_pairs.size() << "   " << denBar.size() << endl;
	for (int iPair = 0; iPair < (int)sorted_list_of_pairs.size()-1; ++iPair)
		cout << "DENSITY: " << sorted_list_of_pairs.at(iPair).first
				<< " <= q <= " << sorted_list_of_pairs.at(iPair+1).first
				<< ": den(q) = " << denBar.at(iPair) << endl;
	*/

		

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
		//cout << "CHECK LHS: " << thisPair.first << "   " << LHS_integral << endl;
	}

	return;
}

double shifter::evaluate_LHS(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			double qz )
{
	int    pairIndex    = 1;
	double previous_qz  = 0.0;
	double this_qz      = 0.0;
	double LHS_integral = 0.0;

	if ( qz < 1.e-20 ) return (0.0);

	while ( this_qz < qz )
	{
		//cout << "checkpoint: " << previous_qz << "   " << this_qz << endl;
		previous_qz   = sorted_list_of_pairs.at(pairIndex-1).first;
		this_qz       = sorted_list_of_pairs.at(pairIndex).first;
		LHS_integral += 2.0 * ( this_qz - previous_qz ) * denBar.at(pairIndex-1);
		pairIndex++;
	}

	//if ( qz > previous_qz )
	//	LHS_integral     += 2.0 * ( qz - previous_qz ) * denBar.at(pairIndex-1);

	return (LHS_integral);
}


//-------------------------------------
// Set righthand side of shift relation.
void shifter::set_RHS(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			vector< pair< double, double > > & RHS )
{

	// RHS similar to LHS
	const int npairs = sorted_list_of_pairs.size();
	//const double one_by_npairs = 1.0 / npairs;
	RHS.reserve( npairs );

	/*for (const auto & eachPair : sorted_list_of_pairs)
	{
		// assumes pair density is constant for now
		const double this_qz = eachPair.first;
		const int this1 = eachPair.second.first;
		const int this2 = eachPair.second.second;

		double RHS_integral = 2.0 * this_qz;

		int npairs_in_average = 0;
		double RHS_BE_enhancement = 0.0;
		for (const auto & iPair : sorted_list_of_pairs)
		{
			const int i1 = iPair.second.first;
			const int i2 = iPair.second.second;

			if ( this1 != i1 and this2 != i2 ) continue;

			Vec4 xDiffPRF = ( allParticles.at(i1).x - allParticles.at(i2).x ) / HBARC;
			//xDiffPRF.bstback( 0.5*(allParticles.at(i1).p + allParticles.at(i2).p) );
			const double Delta_z = xDiffPRF.pz();

			RHS_BE_enhancement += 2.0 * sin(this_qz*Delta_z) / Delta_z;
			++npairs_in_average;
		}

		RHS_integral += RHS_BE_enhancement / static_cast<double>(npairs_in_average);
		RHS.push_back( std::make_pair( eachPair.first, RHS_integral ) );
	}*/

	for (const auto & thisPair : sorted_list_of_pairs)
	{
		const double RHS_integral = evaluate_RHS( sorted_list_of_pairs, RHS, thisPair );
		RHS.push_back( std::make_pair( thisPair.first, RHS_integral ) );
		//cout << "CHECK RHS: " << thisPair.first << "   " << RHS_integral << endl;
	}



	return;
}

double shifter::evaluate_RHS(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			const vector< pair< double, double > > & RHS,
			const pair< double, pair <int,int> > & thisPair )
{
	int    pairIndex = 1;
	double lower_qz  = sorted_list_of_pairs.at(pairIndex-1).first;
	double upper_qz  = sorted_list_of_pairs.at(pairIndex).first;

	const double qz = thisPair.first;

	if ( qz < 1.e-20 ) return (0.0);

	const int this1 = thisPair.second.first;
	const int this2 = thisPair.second.second;

	//cout << sorted_list_of_pairs.size() << "   " << RHS.size() << "   " << pairIndex << endl;
	//cout << "checkpoint(1): " << lower_qz << "   " << upper_qz << "   " << qz << endl;

	// recycle previously computed RHS to save time
	//int dumbcount = 0;
	while ( upper_qz < qz )
	{
		//pairIndex++;
		//cout << "checkpoint(2): " << lower_qz << "   " << upper_qz << "   " << qz << "   " << pairIndex << endl;
		lower_qz = sorted_list_of_pairs.at(pairIndex).first;
		upper_qz = sorted_list_of_pairs.at(++pairIndex).first;
		//dumbcount++;
	}
	//if (dumbcount > 0) pairIndex--;

	//cout << sorted_list_of_pairs.size() << "   " << RHS.size() << "   " << denBar.size() << "   " << pairIndex << endl;
	lower_qz = sorted_list_of_pairs.at(pairIndex-1).first;
	upper_qz = sorted_list_of_pairs.at(pairIndex).first;
	double RHS_integral = RHS.at(pairIndex-1).second;
	//cerr << "Integrating from " << lower_qz << " to " << upper_qz << "; RHS(lower) = " << RHS_integral << endl;

	// the constant piece
	RHS_integral += 2.0 * ( upper_qz - lower_qz ) * denBar.at(pairIndex-1);

	int npairs_in_average = 0;
	double RHS_BE_enhancement = 0.0;

	// the BE enhancement piece
	for (const auto & iPair : sorted_list_of_pairs)
	{
		const int i1 = iPair.second.first;
		const int i2 = iPair.second.second;

		if ( i1<0 or i2<0 ) continue;
		//if ( this1 != i1 and this2 != i2 ) continue;
		if ( this1 != i1 or this2 != i2 ) continue;

		Vec4 xDiff = ( allParticles.at(i1).x - allParticles.at(i2).x ) / HBARC;
		const double Delta_z = xDiff.pz();

		RHS_BE_enhancement += 2.0 * ( sin(upper_qz*Delta_z) - sin(lower_qz*Delta_z) )
							* denBar.at(pairIndex-1) / Delta_z;

		npairs_in_average++;
	}

	RHS_BE_enhancement /= npairs_in_average;
	RHS_integral       += RHS_BE_enhancement;
	pairIndex++;

	//cout << "Result = " << RHS_integral << endl;

	return (RHS_integral);
}


//End of file
