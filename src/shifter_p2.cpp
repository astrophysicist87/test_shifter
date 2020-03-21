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
		evaluate_shift_relation_at_pair( pairs_sorted_by_abs_qz, LHS, RHS, RHS_derivatives );

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
		// ==> THIS IS ALREADY WHAT IS DONE BY DEFAULT

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
	
	void set_LHS_mode2()
	{
		
	}
	
	void set_RHS_mode2()
	{
		
	}



}	//End of namespace

//End of file
