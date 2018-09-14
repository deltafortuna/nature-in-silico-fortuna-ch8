#ifndef POPULATION_H
#define POPULATION_H

#include "params.h"  // provides access to global parameter values (extern)
#include "summarystats.h" // provides access to  summary statistic calculations/functions
/*
using std::find;
using std::sort;
using std::binary_search;
using std::random_shuffle;*/

class Population {

private:
double mu_sequence;
double r_sequence; //////////////////////NEW chapter 6 ////////
vector<Individual*> individuals; 
uniform_int_distribution<int> randompos;
uniform_int_distribution<int> randomind;
uniform_real_distribution<double> randomnum;
poisson_distribution<int> randommut;
poisson_distribution<int> randomrec;                ////// NEW Chapter 6 ////////////
map<int, Allele*> alleles; /// the integer is the position of the allele
ofstream allele_file;
ofstream execution_file;
ofstream sumstat_file;
ofstream abf;    /// new ch 7.4 /// 
ofstream sinfo;  // ch 8.4.1 ////
string rep;
int popn;   //// NEW Chapter 7 /////
bool extant;  //// NEW Chap 7 //
bool activeselection{}, standingvar{}, newvariant{}, nfdselection{}, negselection{}; // set to false initially 
int keypos = -999; // NEW ch 8.4.1
vector<double> fitness; // New ch 8.4.1
double nfd_s, nfd_h; // new ch 8.4.1
vector<int> purifying_sites; // ch 8.4.2

chrono::high_resolution_clock::time_point t1;
chrono::high_resolution_clock::time_point t2;
chrono::high_resolution_clock::time_point r1;
chrono::high_resolution_clock::time_point r2;
chrono::high_resolution_clock::time_point t3;
chrono::high_resolution_clock::time_point t4;
chrono::high_resolution_clock::time_point t5;
chrono::high_resolution_clock::time_point t6;
chrono::high_resolution_clock::time_point t7;
chrono::high_resolution_clock::time_point t8;
chrono::high_resolution_clock::time_point t9;

///// new function 8.4.1
void update_selected_freqAndFit(const int &gen) {
	double count = pop_schedule[popn][gen]*2; // start count at ALL alleles derived 
	vector<double> genotype_freqs = {0.,0.,0.}; // holds frequency of A/A, A/a, a/a, where A is derived allele
	for (auto iter = individuals.begin(); iter != individuals.end(); ++iter) {
		int num_ancestral_alleles = (**iter).get_genotype(keypos);
		count -= num_ancestral_alleles;
		genotype_freqs[num_ancestral_alleles]++;
	}
	double p = count/(pop_schedule[popn][gen]*2); // frequency of derived allele
	double q = 1-p;
	for (int i=0; i<3; ++i) 
		genotype_freqs[i] /= pop_schedule[popn][gen]; 
	cout << gen << ": " << p << endl;
	sinfo << gen << "\t" << p << "\t" << genotype_freqs[0] << "\t" << genotype_freqs[1] << "\t" << genotype_freqs[2] << "\t";

	if (nfdselection) {
		/*
		fitness[0] = 1 - (nfd_s*p*p);
		fitness[1] = 1 - (2*nfd_s*nfd_h*p*q);
		fitness[2] = 1 - (nfd_s*q*q);*/

		double max_fit = -999; 
		fitness[0] = 1 - nfd_h*genotype_freqs[1] + nfd_h*genotype_freqs[2]; 
		if (fitness[0] > max_fit) max_fit = fitness[0];
		fitness[1] = 1 - nfd_s*genotype_freqs[1];
		if (fitness[1] > max_fit) max_fit = fitness[1];
		fitness[2] = 1 - nfd_h*genotype_freqs[1] + nfd_h*genotype_freqs[0];
		if (fitness[2] > max_fit) max_fit = fitness[2];
		for (int i=0; i<3; ++i) 
			fitness[i] /= max_fit;
	}
	sinfo << fitness[0] << "\t" << fitness[1] << "\t" << fitness[2] << endl;
}
//// END new function 8.4.1

bool update_alleles(const int &gen) {  // ch 8.4.1 changed void to bool 
	bool fixtest = true; // ch 8.4.1
	// reset all counts to zero
	for (auto iter = alleles.begin(); iter != alleles.end(); ++iter)
		(*(iter->second)).set_count(0);
	map<int, int> new_allele_counts;

	for (auto iter = individuals.begin(); iter != individuals.end(); ++iter) { 	
		for (int i=0; i<2; ++i) { 

			const vector<int> &a =  (**iter).get_seq(i);  // (1) most recent
			for (auto iter2 = a.begin(); iter2 != a.end(); ++iter2) { //(2) most recent /// THis IS costing the time 			
				++new_allele_counts[*iter2];	// (3) most recent
	//					(*alleles[*iter2]).increment_count();  // NOTE: alleles[*iter] returns a reference (mapped_type) to the pointer to the Allele object at position *iter
								// leading * dereferences the pointer, giving us access
			}			
		}			
	}

	for (auto iter3 = new_allele_counts.begin(); iter3 != new_allele_counts.end(); ++iter3) 
		(*alleles[iter3->first]).set_count(iter3->second);   
	

	// identify lost and fixed alleles and print to allele history file
	vector<int> to_remove;
	for (auto iter = alleles.begin(); iter != alleles.end(); ++iter) {
		int current_count = (*(iter->second)).get_count();

		if (current_count == 0) { // allele LOST from population			
			to_remove.push_back(iter->first);  // first is position 
			if ( iter->first == keypos) fixtest = false; // selective target was lost, restart // NEW ch 8.4.1
		/*	int birthgen = (*(iter->second)).get_birthgen();
			allele_file << iter->first << " " << birthgen << " " << gen - birthgen << " 0" << endl; */
		}
		// all of this BEGIN NEW ch 8.4.1
		if (current_count == pop_schedule[popn][gen]*2) {  // derived allele FIXED in population      // pop_schedule updated in chp 7	
			if ( (activeselection && iter->first != keypos) || !activeselection) { /// New ch 8.4.1   Keep at 1 or everything selected against, may need to change if
				to_remove.push_back(iter->first);  			
				for (auto iter2 = individuals.begin(); iter2 != individuals.end(); ++iter2)
					(**iter2).remove_fixed_allele(iter->first); // removed fixed allele's position from all individuals' sequences (currently stored in nextgen)
			}
		}   /// END NEW largely different at ch 8.4.1 
	}

	// free associated with lost/fixed alleles and remove entry from alleles container
	for (auto iter = to_remove.begin(); iter != to_remove.end(); ++iter) {
		delete alleles[*iter];  // free memory from Allele object itself
		alleles.erase(*iter);  // erase alleles map entry corresponding to the deleted allele
	}		
	
	return(fixtest); // new ch 8.4.1
}

void get_sample(int gen) {

	// new ch 7.3
	ofstream sequencefile;  // will only be created and used if gen % printhapfreq == 0
	string ofname = "deme" + to_string(popn) + "_" + to_string(gen);
	bool printhap=false;
	//if (diploid_sample && gen % printhapfreq == 0) printhap = true; /// oLD conditional only allowed prinhap if diploid_sampling
	if (gen % printhapfreq == 0) printhap = true;    
	if (printhap) sequencefile.open(ofname.c_str());
	/// new ch 7.3 
	vector<bitset<bitlength>> sample; 
	map<int, int> allele_counts; // note that a map is always sorted by keys
	int count = 0;

// new ch 7.3
	int additional = sampsize;
	if (diploid_sample)
		additional /= 2;
// end ch 7.3
cout << "size of individuals vector: " << individuals.size() << endl;
	for (auto iter = individuals.begin(); iter != individuals.begin()+additional; ++iter) { // ch 7.3 additional instead of sampsize
		vector<int> haplotype = (**iter).get_sequence(0);
		for (auto iter2 = haplotype.begin(); iter2 != haplotype.end(); ++iter2) 
			++allele_counts[*iter2];
		/// new ch 7.3
		if (diploid_sample) {
			vector<int> haplotype = (**iter).get_sequence(1);
			for (auto iter2 = haplotype.begin(); iter2 != haplotype.end(); ++iter2) 
				++allele_counts[*iter2];
		}
		// end ch 7.3
	}	

	vector<int> positions;
	for (auto iter = allele_counts.begin(); iter != allele_counts.end(); ++iter) 
 		positions.push_back(iter->first);

 	//// new ch 7.3 ///
 	// print column headers
 	if (printhap) {
 		for (auto iter = positions.begin(); iter != positions.end(); ++iter)
 			sequencefile << "nt" << to_string(*iter) << " ";  
 		sequencefile << endl;
/// new 7.4
 		for (auto iter = positions.begin(); iter != positions.end(); ++iter)
 			sequencefile << (*(alleles[*iter])).get_originating_population() << " ";
 		sequencefile << endl;
 		for (auto iter = positions.begin(); iter != positions.end(); ++iter)
 			sequencefile << (*(alleles[*iter])).get_birthgen() << " ";
 		sequencefile << endl;
/// end new 7.4
 	}
 	/// end new ch 7.3 /// 

	for (auto iter = individuals.begin(); iter != individuals.begin()+additional; ++iter) {  // ch 7.3 additional instead of sampsize
		for (int h=0; h<diploid_sample+1; ++h) {  // new ch 7.3   ; h just 0 if a haploid sample and 0&1 if a diploid sample
			vector<int> haplotype = (**iter).get_sequence(h);   // modified ch 7.3  
			sort(haplotype.begin(), haplotype.end());
			string hap;
			for (auto iter = allele_counts.begin(); iter != allele_counts.end(); ++iter)
				if ( binary_search (haplotype.begin(), haplotype.end(), iter->first)) {   
					hap += "1";
					if (printhap) sequencefile << "1 ";       
				} else {
					hap += "0";
					if (printhap) sequencefile << "0 ";
				}
			sample.push_back(bitset<bitlength> (hap));
			if (printhap) sequencefile << endl;
		}   // end modifications for ch 7.3
	}
	int S = allele_counts.size();

	if (getWindowStats) {
		map<string, vector<double> > stats = get_windowStats(positions, sample, S);
		for (auto iter = stats.begin(); iter != stats.end(); ++iter) {
			sumstat_file << gen << " ";
			sumstat_file << iter->first << " ";
			for (auto iter2 = (iter->second).begin(); iter2 != (iter->second).end(); ++iter2)
				sumstat_file << *iter2 << " ";
			//for (auto i : veck) 
			//	sumstat_file << i << " ";
			sumstat_file << endl;
		}	
	} else {
		double pi = get_pi(sample);
		double watterson = get_watterson(S);
		double tajimasd = get_tajimas_d(pi, watterson, S);
		sumstat_file << gen << " " << pi << " " << watterson << " " << tajimasd << endl;
	}

	if (printhap) sequencefile.close();   // new ch 7.3

/*

t1 = high_resolution_clock::now();	
	pi_file << gen << " " << get_pi(sample) << endl;
t2 = high_resolution_clock::now();
auto duration1 = duration_cast<microseconds>(t2-t1).count();
execution_file << duration1 << " ";
t1 = high_resolution_clock::now();	
	watterson_file << gen << " " << get_watterson(sample, S) << endl;
t2 = high_resolution_clock::now();
auto duration2 = duration_cast<microseconds>(t2-t1).count();
	tajimasd_file << gen << " " << get_tajimas_d()
execution_file << duration2 << endl;
*/
}

vector<vector<int> > mutate(const vector<int> &parents, const int &gen) {
	vector<vector<int> > mutation_results;

	// determine which, if any, positions are mutated    ///// NEW
	vector<int> mutnum{randommut(e)};  //// NEW
	mutnum.push_back(randommut(e));  //// NEW

	mutation_results.push_back({mutnum[0]});
	mutation_results.push_back({mutnum[1]});
	
	// determine which of the two homologs is transmitted by each parent
	mutation_results.push_back({(randomnum(e)<0.5) ? 0 : 1}); 
	mutation_results.push_back({(randomnum(e)<0.5) ? 0 : 1});

	// resolve any mutation(s) that did occur
	for (int i=0; i<2; ++i)  {  /// NEW 
		for (int j = 0; j < mutnum[i]; ++j)  { // loop not entered if no mutation (i.e., mutnum[i] == 0)
			int position = randompos(e);
			if (alleles.find(position) == alleles.end()) { // new mutation to a derived allele in the population
				alleles.insert({position, new Allele(position, gen, popn)});
//				if (trackAlleleBirths) 
//					abf << "nt" << position << "\t" << gen << "\t" << popn << endl;   /// new ch 7.4  /// TRUE COMMENT: allele birth file
				mutation_results[i].push_back(position);
			} else { // mutation present in POPULATION; determine if derived allele found in the considered sequence
				vector<int> seq = (*(individuals[parents[i]])).get_sequence(mutation_results[i+2][0]); 
				vector<int>::iterator p = find(seq.begin(), seq.end(), position);
				if (p != seq.end())   // back mutation
					mutation_results[i].push_back(position *  -1);	  // negative position signals removal of allele by back mutation
				else
					mutation_results[i].push_back(position);
			}	
 		}
	}

	return mutation_results;
}

vector<int> recombine() {
	vector<int> breakpoints;
	if (useRec) {  /// if false, empty vector passed to Individual constructor
 		int chiasmata = randomrec(e);
 		for (int i=0; i<chiasmata; ++i) {
	 		if (useHotRec) {
	 			double c = randomnum(e);
	 			if (c < hotrecStart * recrate / r_sequence )
	 				breakpoints.push_back( c * r_sequence / recrate );
	 			else if (c < (hotrecStop*hotrecrate - hotrecStart*(hotrecrate-recrate)) / r_sequence)
	 				breakpoints.push_back(  (c*r_sequence + hotrecStart*(hotrecrate - recrate))    /     hotrecrate);
	 			else
	 				breakpoints.push_back( (c*r_sequence - (hotrecStop-hotrecStart)*(hotrecrate-recrate)   )    / recrate );
		 	} else {
				breakpoints.push_back(randompos(e));
			}
		}
		sort(breakpoints.begin(), breakpoints.end());
	}

	return breakpoints;
}

public:

/// NEW CHapter 7 /// 
inline int get_popnum() { return popn;}   
inline bool get_extant() {return extant;}
inline int get_individuals_vectorsize() {return individuals.size();}
//// new chapter 7.3
inline int get_current_popsize(int gen) {return pop_schedule[popn][gen];}   

vector<int> set_extant() {
	extant = 1;
	vector<int> i;
	if (splitgenesis[popn][0] > 0) {
		i.push_back(1);
		i.push_back(splitgenesis[popn][1]); // source population
		i.push_back(splitgenesis[popn][2]); // percent 
	} else if (mergegenesis[popn][0] > 0 ) {
		i.push_back(2);
		i.push_back(mergegenesis[popn][1]); // first source population
		i.push_back(mergegenesis[popn][2]); // second source population
	} else
		i.push_back(0);
	return(i);
}
//// new chp 7.3

inline void set_extinct() {extant = 0;}
inline vector<vector<int>> get_sequences(int indnum) { return (*individuals[indnum]).get_sequences();}
inline void add_immigrant(vector<vector<int>> ses) {individuals.push_back( new Individual(ses) ); }

void remove_emigrants(int Nm) {
	for (auto iter = individuals.begin(); iter != individuals.begin() + Nm; ++iter)
		delete *iter;
	individuals.erase(individuals.begin(), individuals.begin()+Nm);
}

vector<int> get_allele_positions() {
	vector<int> v;
	for (auto iter=alleles.begin(); iter!=alleles.end(); ++iter)
		v.push_back(iter->first);
	return v;
}

vector<int> get_allele_info(int s) {
	vector<int> v = {s};
	v.push_back( (*alleles[s]).get_birthgen() );
	v.push_back( (*alleles[s]).get_originating_population() );
	return v; 
}

void insert_new_allele(vector<int> v) {
	alleles.insert( { v[0]  , new Allele( v[0], v[1], v[2] ) } );
}
//// END NEW Chapter 7 ////////

void reproduce(int gen) {
	// reparameterize uniform distribution based on current population size	
	randomind.param(uniform_int_distribution<int>::param_type(0, pop_schedule[popn][gen==0 ? gen : gen-1] - 1) );     // pop_schedule updated in chp 7  // early7 WAS IT? WHY? CHANGED BACK; HAD TO ADD THE -1 because chance you get highest number on very first roll, which is out of bounds
//cout << "migration detail: " << individuals.size() << " is size of individuals vector for population #" << popn << " before anything." << endl;		
	/// early7 ///
	int N = pop_schedule[popn][gen]; // early7 // how many individuals to make for next gen, taking into consideration number of immigrants and emigrants
	if (modelMigration) {
		int n_imm = 0;
		int n_emi = 0; 
		for (int i=0; i<pop_num; ++i) {
//cout << "popn #" << i << endl;			
			n_emi += mig[popn][i] * pop_schedule[i][gen];
			n_imm += mig[i][popn] * N;
		}
		N += n_emi;
		N -= n_imm;
//cout << "migration detail: " << N << " is calculated N for population #" << popn << endl;		
	}
	///// end early7 /////

	///NEW ch 8.4.1 ///																	// should be for (int i=0; i< pop_schedule[popn][gen]; ++i), because this is the number to make for current generation
	// removal of potential parents by selection
	vector<int> recode_parents;
	if (activeselection) 
		for (int i=0; i<individuals.size(); ++i) 
			if (randomnum(e) <= fitness[ (*individuals[i]).get_genotype(keypos) ] ) recode_parents.push_back(i); 
	if (negselection) {
vector<int> geno_counts = {0, 0, 0};
		for (int i=0; i<individuals.size(); ++i) {
			if (randomnum(e) <= fitness[ (*individuals[i]).get_negsel_genotype(purifying_sites) ] ) recode_parents.push_back(i);
//geno_counts[ (*individuals[i]).get_negsel_genotype(purifying_sites) ]++;
		}
//cout << "At generation " << gen << ", genotype counts: " << geno_counts[0] << " " << geno_counts[1] << " " << geno_counts[2] << endl;
	}	
	
	if (activeselection || negselection)	{
		randomind.param(uniform_int_distribution<int>::param_type(0, recode_parents.size() - 1) ); 
//cout << "recode size is: " << recode_parents.size() << endl; 
	}
	/// End NEW ch 8.4.1

	#pragma omp parallel for
	for (int i=0; i<N; ++i) {  /// early 7
		vector<int> parents;
		///NEw ch 8.4.1
		if (activeselection  || negselection) {//selection
			parents.push_back(recode_parents[randomind(e)]);
			parents.push_back(recode_parents[randomind(e)]);
		} else { // end new ch 8.4.1
			parents.push_back(randomind(e));
			parents.push_back(randomind(e));
//cout << i << ": " << parents[0] << " " << parents[1] << endl;			
		}
		// create descendant of individuals parents[0] and parents[1]
		individuals.push_back( new Individual(individuals[parents[0]], individuals[parents[1]], mutate(parents, gen), recombine() ) );
	}
	// delete dynamically allocated individuasl of the last generation
//cout << "migration detail: " << individuals.size() << " is size of indiv vector after repro but before erasure for pop #" << popn << endl;
	for (auto iter = individuals.begin(); iter != individuals.end() - N; ++iter) //- pop_schedule[popn][gen==0 ? gen : gen-1]; ++iter)  // pop_schedule updated in chp 7
		delete *iter;
	individuals.erase(individuals.begin(), individuals.end() - N);  //  - pop_schedule[popn][gen==0 ? gen : gen-1]);  // remove the orphaned pointers from individuals   // pop_schedule updated in chp 7
//cout << "migration detail: " << individuals.size() << " is size of indiv vector after reproduction AND erasure for pop #" << popn << endl;

	if (activeselection) update_selected_freqAndFit(gen);   // new 8.4.1
}

//// early7 /// Made its own function, because need to wait until after migration to perform sampling !!!!! 
bool sample(int gen) {  // ch 8.4.1 changed void to bool 
	bool fixtest = update_alleles(gen);  // ch8.4.1 added bool
	if (!fixtest && sellocus[popn][4] == 1) {return fixtest;}  // ch 8.4.1, worthwhile skipping sample calcualtion if selected target lost
	random_shuffle(individuals.begin(), individuals.end() ) ;
	get_sample(gen +1 ); //// IF this works change in earlier exposition, not in ch 7.3 (though this is where i thought of it)		
	return(fixtest);  // new ch 8.4.1
}	
//// early7  /// 

void close_output_files () {
	//allele_file.close(); 
	sumstat_file.close();
	execution_file.close();
	if (trackAlleleBirths) abf.close();   /// abf ////  addded to text of ch 7.3 
	if (activeselection) sinfo.close(); // ch 8.4.1
}

Population (int popnum, int eextant):popn(popnum), extant(eextant) {  //// NEW Chapter 7, specifically the artgument popnum

	rep  = "ch8pt4-pt1";	

	randompos.param(uniform_int_distribution<int>::param_type(1,seqlength));
	randomind.param(uniform_int_distribution<int>::param_type(0,pop_schedule[popn][0] - 1));   // pop_schedule updated in chp 7  ///early7 does this need the -1 ???
	randomnum.param(uniform_real_distribution<double>::param_type(0.,1.));
	mu_sequence = seqlength * mutrate; 
	randommut.param(poisson_distribution<int>::param_type(mu_sequence)); 

	if (useRec) {
		if (useHotRec) {
			int hotspot_length = hotrecStop - hotrecStart + 1;
			r_sequence = ( hotspot_length * hotrecrate )  +  ((seqlength - hotspot_length) * recrate);
		} else 
			r_sequence = seqlength * recrate;
		randomrec.param(poisson_distribution<int>::param_type(r_sequence));     
	}

     // New ch 8.4.1
	if (possel[popn][0] != 0 || nfdsel[popn][0] !=0 || negsel[popn][0] != 0) {
		if (nfdsel[popn][0] != 0) 
			nfdselection = true;
		if (negsel[popn][0] != 0) {
			negselection = true;
			fitness = {1, 1-(negsel[popn][0]*negsel[popn][1]), 1- negsel[popn][0]};
			vector<int> possible_sites(negsel[popn][3]-negsel[popn][2]+1);
			iota(possible_sites.begin(), possible_sites.end(), negsel[popn][2]);
			random_shuffle(possible_sites.begin(), possible_sites.end());
			auto iter = possible_sites.begin();
			purifying_sites.assign(iter, iter+negsel[popn][4]); 
			sort (purifying_sites.begin(), purifying_sites.end());
for (auto zz = purifying_sites.begin(); zz != purifying_sites.end(); ++zz) {
	cout << *zz << endl;
}				
		} else 
			activeselection	= true;

		if (activeselection) {
			keypos = sellocus[popn][0];  // temporary if selection on standing variation
			if (sellocus[popn][1] > 1) 
				standingvar=true;
			else 
				newvariant=true;
			double s = possel[popn][0];
			double t = possel[popn][1];
			double h = possel[popn][2]; 
			nfd_s = nfdsel[popn][0];
			nfd_h = nfdsel[popn][1];
			if (t == 0) fitness = {1, 1-(h*s), 1-s}; // dominance or additivity 
			else fitness = {1-s, 1, 1-t};		 // overdominance
		}
	} 
	// End New ch 8.4.1
cout << "standingvar is: " << standingvar << endl;
cout << "activeselection is " << activeselection << endl;
cout << "nfdselection is: " << nfdselection << endl;
cout << "keypos is: " << keypos << endl;
cout << "negselection is: " << negselection << endl;
	individuals.reserve(popsize*10);

	// ch 7.4 // 
	string ofname = "deme" + to_string(popn) + "_allele_births";
	if (trackAlleleBirths)
		abf.open(ofname.c_str());
	/// ch 7.4 ///   

	/// New ch 8.4.1
	if (activeselection) {
		ofname = "deme" + to_string(popn) + "_selectioninfo";
		sinfo.open(ofname.c_str());
		sinfo << "Standing Variation? " << standingvar << endl;
		sinfo << "Negative frequency-dependence? " << nfdselection << endl;
	}
	/// END new ch 8.4.1

	if (useMS[popn]) { // start population with MS generated variation  // index to useMS added chp 7
		cout << "using MS to initialize population ..." << endl;
		system(mscommand[popn].c_str());   /// added [popn] in Chp 7

		ifstream ms_output("ms_output");
		string ms_line;
		regex query("positions");
		bool trigger = false;
		vector<int> allele_positions;

	 	//// NEW ch 8.4.1
		if (activeselection && !standingvar)     ///// NECESSARY ?????????????
 			//if (! (find(allele_positions.begin(), allele_positions.end(), keypos) != allele_positions.end()) ) { 
	 		//allele_positions.push_back(keypos);  // add position of selected locus
	 		alleles.insert( { keypos, new Allele(keypos, -1, popn) } );
	 		//}	
 			//// END NEW ch 8.4.1	
		
	
 		while(getline(ms_output, ms_line)) {   //use until loop first??
 												// should obviate need for trigger
			if (regex_search(ms_line, query)) {
				trigger = true;
				cout << "got it" << endl;
				istringstream iss(ms_line);
				string s;
				iss >> s; //skip the first subpart, which is "positions:"
				while (iss >> s) { // read decimal positions, convert to base pair position,
									// and create new allele at that positon
					int position = seqlength * atof(s.c_str());
					allele_positions.push_back(position);
				//	if (!standingvar && position == keypos )  /// new ch 8.4.1
				//		continue;   // new ch 8.4.1
				//	else { // new ch 8.4.1
					if (standingvar || position != keypos)
					//	allele_positions.push_back(position);
						alleles.insert( { position  , new Allele(position,-1, popn) } ); // the origin generation unknown, (dummy -1)
					//}  // new ch 8.4.1
				}						
				continue;
 			}	

			if (trigger) { // then allele positions have been determined and each line is a haplotype
 				vector<int> s1, s2;
 				for (int i=0; i < ms_line.length(); ++i) 
 					if (ms_line[i] == '1') {
 						if (allele_positions[i] != keypos) // NEW ch 8.4.1 only ignore if MS happened to include variant at target of selection
 							s1.push_back(allele_positions[i]);  // ch 8.4.1
 						if (standingvar) 
 							(*(alleles[allele_positions[i]])).increment_count();   // ch 8.4.1
 					}		
 				getline(ms_output, ms_line);
 				for (int i=0; i < ms_line.length(); ++i) 
 					if (ms_line[i] == '1') {
 						if (allele_positions[i] != keypos) // NEW ch 8.4.1
 							s2.push_back(allele_positions[i]);  // ch 8.4.1
 						if (standingvar) 
 							(*(alleles[allele_positions[i]])).increment_count();  // ch 8.4.1
 					}	
 				/// New ch 8.4.1
 				if (activeselection && !standingvar  && newvariant) { 
 					s1.push_back(keypos);
 					newvariant = false;
 					sort(s1.begin(), s1.end());
 				}
 				//if (activeselection && !standingvar) {
	 			//	sort(s1.begin(), s1.end()); // Move selected position from back of sequence to its proper place
	 			//	sort(s2.begin(), s2.end());  //not necessary anymore; AND can be moved up to previous loop
 				//}  				

 				/// End New ch 8.4.1
 				vector<vector<int>> ses{s1,s2};
 				individuals.push_back( new Individual(ses) );
 			}
		}
	} else {
		for (int NN=0; NN<10000; ++NN) {   // have to fix!!!! so can use with other stuff
			vector<int> s1; vector<int> s2;
			vector<vector<int> > ses{s1,s2};
			individuals.push_back( new Individual(ses) );
		}

	}
//// NEw 8.4.1
	if (activeselection && standingvar) {   // then need to determine which locus roughly meets criteria to become keypos
		int low=keypos-floor(sellocus[popn][2]*seqlength);
		if (low < 0 ) low = 0;
		int high=keypos+floor(sellocus[popn][2]*seqlength);
		if (high >= seqlength) high = seqlength-1;
		vector<int> current_best{-999, 999999,-999};
		int acceptablediff = floor(pop_schedule[popn][0]*2*sellocus[popn][3]); 
		int target_low = sellocus[popn][1] - acceptablediff;
		int target_high = sellocus[popn][1] + acceptablediff;	
		for (auto iter=alleles.begin(); iter != alleles.end(); ++iter) {
			if (iter->first >= low && iter->first <= high) { // in range
				int ct = (*(iter->second)).get_count(); 			
				if (ct >= target_low && ct <=target_high) {
					int diff = abs(sellocus[popn][1]-ct);				
					if (diff < current_best[1] ) { 						
						current_best[0] = iter->first;
						current_best[1] = diff;
						current_best[2] = ct;						
					}
				}
			} else {
				if (iter->first > high) break;
				else continue;
			}
		}
		if (current_best[0] != -999) {
			keypos = current_best[0];	
			sinfo << "Standing variant identified at base pair " << keypos << " with a starting allele count at time 0 of " << current_best[2] << endl;
		} else {
			cout << "Did not find a suitable standing variant." << endl;
			throw("suitable standing variant not identified");
		}	

		if (nfdselection) 
 			update_selected_freqAndFit(0);

 		sinfo << "gen\tp\tf_AA\tf_Aa\tf_aa\tw_AA\tw_Aa\tw_aa" << endl;
	}
//// END 8.4.1
cout << "keypos is " << keypos << endl;
	string fname{"allele_history"};

	fname = "sumstats" + rep + "." + to_string(popn);
	sumstat_file.open(fname.c_str());
///// new CH6 ////
	if (getWindowStats) {
		sumstat_file << "gen stat ";
		for (int w=0; w + windowSize <= seqlength; w += windowStep) 
			sumstat_file << "w" << w+windowStep << " "; // window name labels the ending position of the window
		sumstat_file << endl;
	}	else 	
		sumstat_file << "gen pi watterson tajimasd" << endl; 
/////// new ch6 ///////////////

	/*fname = "execution" + rep;
	execution_file.open(fname.c_str());
	execution_file << "pi.time wat.time" << endl;
	*/
}

static mt19937 e;

};

#endif
