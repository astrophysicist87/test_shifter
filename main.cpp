#include <iostream>

#include "include/shifter.h"
#include "include/ParameterReader.h"
#include "include/random_events.h"

using namespace std;
using namespace shift_lib;

int main(int argc, char *argv[])
{
	// Read-in free parameters
	ParameterReader * paraRdr = new ParameterReader;
	paraRdr->readFromFile("./parameters.dat");

	// Read-in command-line arguments
	paraRdr->readFromArguments(argc, argv);
	paraRdr->echo();

	// Vector to hold all event information
	vector<ParticleRecord> allParticles;

	// Read in the files
  cout << "Generating events...";
	generate_events(allParticles, paraRdr);
  cout << "finished." << endl;

	// Create HBT_event_generator object from allEvents
	shifter event( paraRdr, allParticles, cout, cerr );

	cerr << "Finished 1 event" << endl;

	///*
	// Loop over several events
	const int nLoops = paraRdr->getVal("RNG_nLoops");
	for (int iLoop = 1; iLoop < nLoops; ++iLoop)
	{
		generate_events(allParticles, paraRdr);
		event.update_records( allParticles );
		cerr << "Finished "1 << iLoop+1 << " events" << endl;
	}
	//*/

	//for ( auto const particle: allParticles )
	//	cout << particle;

	// Done.
	return (0);
}

//End of file
