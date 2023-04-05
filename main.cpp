#include <iostream>
#include <vector>

#include <omp.h>

#include "include/shifter.h"
#include "include/ParameterReader.h"
#include "include/random_events.h"

using namespace std;
using namespace shift_lib;

double standard_deviation( const vector<double> & v )
{
	const double n = v.size();
	double x2 = 0.0, runsum = 0.0, variance = 0.0;
	for (const auto e: v)
	{
		variance -= e*runsum;
		runsum += e;
		x2 += e*e;
	}
	variance *= 2./(n*(n-1));
	variance += x2/n;
	return sqrt(variance);
}

int main(int argc, char *argv[])
{
	// Read-in free parameters
	ParameterReader * paraRdr = new ParameterReader;
	paraRdr->readFromFile("./parameters.dat");

	string results_directory = std::string(argv[1]);
	string shift_mode = std::string(argv[2]);

	// Read-in command-line arguments (skip one to require results directory path)
	paraRdr->readFromArguments(argc, argv, (string)("#"), 3);
	paraRdr->echo();

	// Loop over several events
	const int nLoops = paraRdr->getVal("RNG_nLoops");
	int number_of_completed_events = 0;

	#pragma omp parallel for schedule(static)
	for (int iLoop = 0; iLoop < nLoops; ++iLoop)
	{
		// Vector to hold all event information
		vector<ParticleRecord> allParticles;

		// Read in the files
		generate_events(allParticles, paraRdr);

		// Create shifter object for each event
		shifter event( paraRdr, allParticles, shift_mode, cout, cerr );

		#pragma omp critical
		{
			event.print( number_of_completed_events, allParticles, results_directory + "/events.dat" );
			cerr << "Finished " << number_of_completed_events++ << " events" << endl;
		}
	}

	// Done.
	return (0);
}

//End of file
