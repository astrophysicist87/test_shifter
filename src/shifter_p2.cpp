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
	//--------------------------------------------------------------------------
	// Perform Bose-Einstein corrections on an event.

	void shifter::shiftEvent_v2()
	{

		// Loop through pairs of identical particles and find shifts.
		bool enoughPairsToProceed = setSortedPairs( allParticles );
		if ( enoughPairsToProceed )
		{
			cout << "shifterCheck: NPair = " << sortedPairs.size() << endl;
			shiftPairs_mode2();
		}

		// Add in shifts without compensations
		allParticles_Shifted = allParticles;
		for (auto & thisParticle: allParticles_Shifted)
			thisParticle.p += thisParticle.pShift;

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
	


		constexpr bool perform_compensation = false;
		if ( not perform_compensation )
			cout << "WARNING: compensation currently turned off!" << endl;


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


		constexpr bool check_for_bad_events = false;
		if ( not check_for_bad_events )
			cout << "WARNING: checking for bad events currently turned off!" << endl;


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

		// Done.
		return;

	}
	
	void shifter::shiftPairs_mode2()
	{
		vector< pair< double, double > > LHS, RHS, RHS_derivatives;
		evaluate_shift_relation_at_pair_mode2( pairs_sorted_by_abs_qz, LHS, RHS, RHS_derivatives );

		compute_shifts( pairs_sorted_by_abs_qz, LHS, RHS, RHS_derivatives );

		return;
	}
	
	

	//-------------------------------------
	// Compute the unshifted pair integrals
	// at each pair's relative momentum.
	void shifter::evaluate_shift_relation_at_pair_mode2(
				const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
				vector< pair< double, double > > & LHS,
				vector< pair< double, double > > & RHS,
				vector< pair< double, double > > & RHS_derivatives
				)
	{
		// -----------------------
		// Shift relation - mode 2
		// -----------------------
		// n - number of particles
		// N = n(n-1)/2 - number of pairs
		// Instead of shifting N pairs and using to compute n particle shifts (overdetermined),
		// compute ONLY n shifts, one for each particle, computed from the average of all the
		// pairs to which it belongs
		// --> but shifts are only computed for a given pair, not for individual particles!
		// --> how to directly compute shift for only single particle?
		// ==> BEST WAY IS APPARENTLY JUST TO GET N SHIFTS AND COMPUTE AVERAGE SHIFT FOR EACH PARTICLE
		// ==> THIS IS ALREADY WHAT IS DONE BY DEFAULT...

		//-------
		// Reset.
		LHS.clear();
		RHS.clear();
		RHS_derivatives.clear();

		// Make sure pair density is stored, if needed.
		set_pair_density( sorted_list_of_pairs );

		//------------------
		// Set LHS integral.
		set_LHS_mode2( sorted_list_of_pairs, LHS );

		//------------------
		// Set RHS integral.
		set_RHS_mode2( sorted_list_of_pairs, RHS, RHS_derivatives );

		return;
	}
	
	void shifter::set_LHS_mode2(
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
	
	void shifter::set_RHS_mode2(
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
			const double RHS_integral = evaluate_RHS_mode2( sorted_list_of_pairs, thisPair, thisPair.first, RHS_derivative );
			RHS.push_back( std::make_pair( thisPair.first, RHS_integral ) );
			RHS_derivatives.push_back( std::make_pair( thisPair.first, RHS_derivative ) );
			cout << "CHECK RHS: " << thisPair.first << "   " << RHS_integral << endl;
		}

		return;
	}


	double shifter::evaluate_RHS_mode2(
				const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
				const pair< double, pair <int,int> > & thisPair,
				const double qz_in, double & RHS_derivative )
	{
		int    pairIndex = 1;
		double lower_qz  = sorted_list_of_pairs.at(pairIndex-1).first;
		double upper_qz  = sorted_list_of_pairs.at(pairIndex).first;

		double qz = abs(qz_in);
		if ( qz < 1e-20 ) return (0.0);

		const int this1 = thisPair.second.first;
		const int this2 = thisPair.second.second;

		double RHS_integral = 0.0;

		while ( upper_qz < qz )
		{
			// the constant piece (integrand assumed to be symmetric; hence factor of 2)
			RHS_integral += 2.0 * ( upper_qz - lower_qz ) * denBar.at(pairIndex-1);
			//RHS_derivative = 2.0 * denBar.at(pairIndex-1);

			int npairs_in_average = 0;
			double RHS_BE_enhancement = 0.0;
			//double RHS_BE_enhancement_derivative = 0.0;
	
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
	
				RHS_BE_enhancement += 2.0 * ( sin(upper_qz*Delta_z) - sin(lower_qz*Delta_z) )
									* denBar.at(pairIndex-1) / Delta_z;
	
				//RHS_BE_enhancement_derivative += 2.0 * cos(upper_qz*Delta_z) * denBar.at(pairIndex-1);
	
				npairs_in_average++;
			}
	
			RHS_BE_enhancement /= npairs_in_average;
			//RHS_BE_enhancement_derivative /= npairs_in_average;
			RHS_integral += RHS_BE_enhancement;
			//RHS_derivative += RHS_BE_enhancement_derivative;

			// shift to next interval
			pairIndex++;
			lower_qz = sorted_list_of_pairs.at(pairIndex-1).first;
			upper_qz = sorted_list_of_pairs.at(pairIndex).first;
		}

		// ------------------------
		// Do last partial interval
		// ------------------------

		// the constant piece (integrand assumed to be symmetric)
		RHS_integral += 2.0 * ( qz - lower_qz ) * denBar.at(pairIndex-1);
		RHS_derivative = 2.0 * denBar.at(pairIndex-1);

		// Use a block here in order to re-use variables used above
		{
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
		}

		return (RHS_integral);
	}


}	//End of namespace

//End of file
