#ifndef SHIFTER_H
#define SHIFTER_H

#include <omp.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <complex>

#include "ParameterReader.h"
#include "ParticleRecord.h"

using namespace std;

class shifter
{
	private:

		ParameterReader * paraRdr;

		vector<ParticleRecord> allParticles;
		
		// miscellaneous
		string path;
		ostream & out;
		ostream & err;

		bool include_pair_density;

		// The pair density
		vector<double> denBar;

		// Various ways to sort pairs
		vector< pair< double, pair <int,int> > > sortedPairs;	// sorted by Q^2
		vector< pair< double, pair <int,int> > > pairs_sorted_by_qzPRF;
		vector< pair< double, pair <int,int> > > pairs_sorted_by_abs_qzPRF;
		vector< pair< double, pair <int,int> > > pairs_sorted_by_qz;
		vector< pair< double, pair <int,int> > > pairs_sorted_by_abs_qz;



	public:

		// Constructors, destructors, and initializers
		shifter( ParameterReader * paraRdr_in,
					const vector<ParticleRecord> & allParticles_in,
					ostream & out_stream = std::cout,
					ostream & err_stream = std::cerr )
					:
					out(out_stream),
					err(err_stream)
					{ initialize_all( paraRdr_in, allParticles_in ); };


		void initialize_all( ParameterReader * paraRdr_in, const vector<ParticleRecord> & allParticles_in );

		void update_records( vector<ParticleRecord> & allParticles_in );

		~shifter();

		void shiftEvent();

		bool setSortedPairs();

		void shiftPairs_mode1();

		void set_LHS(
				const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
				vector< pair< double, double > > & LHS );

		double evaluate_LHS(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			double qz );

		void set_RHS(
				const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
				vector< pair< double, double > > & RHS );

		double evaluate_RHS(
			const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
			const vector< pair< double, double > > & RHS,
			const pair< double, pair <int,int> > & thisPair,
			double & RHS_derivative );

		void compute_shifts( const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs );

		void evaluate_shift_relation_at_pair(
				const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs,
				vector< pair< double, double > > & LHS,
				vector< pair< double, double > > & RHS
				);

		//double Newtons_Method( const double a, const double b );

		void set_pair_density( const vector< pair< double, pair <int,int> > > & sorted_list_of_pairs );

};

// End of file
#endif
