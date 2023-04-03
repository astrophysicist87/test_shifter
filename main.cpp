#include <iostream>
#include <vector>

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

	string results_directory = std::string(argv[1]) + "/";

	// Read-in command-line arguments (skip one to require results directory path)
	paraRdr->readFromArguments(argc, argv, 2);
	paraRdr->echo();

	// Vector to hold all event information
	vector<ParticleRecord> allParticles;

	// Read in the files
	generate_events(allParticles, paraRdr);

	// garbage needly quickly
	vector<double> v0, v1;
	for (auto i: allParticles) v0.push_back(i.p.pz());

	// Create HBT_event_generator object from allEvents
	shifter event( paraRdr, allParticles, cout, cerr );

	for (auto i: allParticles) v1.push_back(i.p.pz());
	event.print( 0, allParticles, results_directory + "events.dat" );


	// Loop over several events
	const int nLoops = paraRdr->getVal("RNG_nLoops");
	for (int iLoop = 1; iLoop < nLoops; ++iLoop)
	{
		generate_events(allParticles, paraRdr);
		for (auto i: allParticles) v0.push_back(i.p.pz());
		event.update_records( allParticles );
		for (auto i: allParticles) v1.push_back(i.p.pz());
		event.print( iLoop, allParticles, results_directory + "events.dat" );
		cerr << "Finished " << iLoop+1 << " events" << endl;
	}

	double s0 = standard_deviation(v0);
	double s1 = standard_deviation(v1);

	//for (int i = 0; i < v0.size(); i++)
	//	cerr << i << "   " << v0.at(i) << "   " << v1.at(i) << "\n";

	cerr << "Unshifted sigma = " << s0 << endl;
	cerr << "Shifted sigma = " << s1 << endl;


	//for ( auto const particle: allParticles )
	//	cout << particle;

	// Done.
	return (0);
}

//End of file
