#include <algorithm>
#include <string>
#include <iterator>
#include <bitset>
#include <ctime>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <memory>
#include <chrono>
#include <regex>
#include <numeric>    // new ch 8.4.2 for iota()

/*
using std::string; // for some strange reason only compiles if using statements placed before the
			// custom header files;  Need to find out why that's the case here!!!!
using std::vector;
using std::map;
using std::list;
using std::cout;
using std::cin;
using std::endl;
using std::time;
using std::ifstream;
using std::ofstream;
using std::bitset;
using std::getline;
using std::istringstream;
using std::mt19937;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::poisson_distribution; 
using std::regex;
using std::regex_search;
//// new chpater 6 ////
using std::upper_bound;
///// new chapter 6 /////

using namespace std::chrono;*/

using namespace std;

#include "matrix.h" // new chapter 7 ////////
#include "params.h"
#include "params.cc"
#include "allele.h"
#include "individual.h"
#include "population.h"
#include "metapopulation.h"  // new chapter 7 /////////////

//int main(int argc, char *argv[]) {
int main() {
	//high_resolution_clock::time_point execution;
	//high_resolution_clock::time_point execution2;

	//int gens = atoi(argv[1]);

	//execution = high_resolution_clock::now();

	mt19937 engine(time(0));  //initialize the random engine
	Population::e  = engine;

	mt19937 engine2(time(0));
	Metapopulation::f = engine2;

	cout << seqlength << " is sequence length" << endl;

	cout << "MIGRATION matrix:" << endl;
	mig.print_matrix();
	bool simulate = true; // new ch 8.4.1
	
	while (simulate) {  /// new ch 8.4.1
		Metapopulation meta;

		for (int i =0; i < runlength; ++i) {
			if (i % 5 == 0) cout << i<< endl;
			bool test = meta.reproduce_and_migrate(i);
			if (! test) break; // new ch 8.4.1
			if (i == runlength - 1) simulate=false; 	 // new ch 8.4.1
		}
		meta.close_output_files();
	}  /// new ch 8.4.1
	//execution2 = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>( execution2 - execution ).count();	

	/*string fname = "execution";
	fname += rep;
	ofstream execution_file;
	execution_file.open(fname.c_str());
	execution_file << duration << endl;
	execution_file.close(); */

	return 0;
}

// static variables for population class 
mt19937 Population::e;
mt19937 Metapopulation::f;