#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "allele.h"

class Individual {

private:
	vector<vector<int> > sequences; 

	void remove_allele_by_position (int seqnum, int position) {
		auto pos = find(sequences[seqnum].begin(), sequences[seqnum].end(), position);
		if (pos != sequences[seqnum].end())  // ensures that allele's position is in the vector
			sequences[seqnum].erase(pos);
	}

	void resolve_crossover (vector<int>& breaks) {
		int numbreaks = breaks.size();
		vector<int>::iterator lower, upper;
		map<int, vector<int> > segments;
		vector<int> newvec;

		for (int seq = 0; seq < 2; ++seq) {
			lower = sequences[seq].begin();
			for (int i=0; i<numbreaks; ++i) {
				upper  = upper_bound(sequences[seq].begin(), sequences[seq].end(), breaks[i]);
				newvec.assign(lower, upper);
				lower = upper;
				segments[i + seq*(numbreaks+1)] = newvec;
			}
			newvec.assign(lower, sequences[seq].end());
			segments[numbreaks + seq*(numbreaks+1)] = newvec;
		}

		for(int i=0; i<= numbreaks; ++i) {
			if (i%2 == 0) {
				if (i == 0) {sequences[0] = segments[0]; sequences[1] = segments[numbreaks+1];}
				else {
					sequences[0].insert(sequences[0].end(), segments[i].begin(), segments[i].end());     
					sequences[1].insert(sequences[1].end(), segments[i+numbreaks+1].begin(), segments[i+numbreaks+1].end());
				}
			} else{
				sequences[0].insert(sequences[0].end(), segments[i+numbreaks+1].begin(), segments[i+numbreaks+1].end());
				sequences[1].insert(sequences[1].end(), segments[i].begin(), segments[i].end());
			}
		}
	}

public:

	inline vector<vector<int> > get_sequences() { return sequences; }

	inline vector<int> get_sequence(int whichseq) { return sequences[whichseq]; } 

	const vector<int> & get_seq(int whichseq) { return sequences[whichseq]; }

	void remove_fixed_allele(int to_remove) {

		for (int i = 0; i<2; ++i) {
			vector<int>::iterator p = find(sequences[i].begin(), sequences[i].end(), to_remove); // assumes only one instance of allele, which is safe in this case
			if (p != sequences[i].end()) // == myVector.end() means the element was not found   
		    	sequences[i].erase(p); 
		 }

	}

/*	vector<int> get_alleles() { 
		vector<int> a = sequences[0];
		a.insert(a.end(), sequences[1].begin(), sequences[1].end());
		return(a);
	}
*/
	//// new ch 8.4.1
	int get_genotype(int pos) { // retturns number of ancestral alleles (i.e., 0s)
		int geno = 2;
		if ( find( sequences[0].begin(), sequences[0].end(), pos) != sequences[0].end() ) --geno;
		if ( find( sequences[1].begin(), sequences[1].end(), pos) != sequences[1].end() ) --geno;
		return geno; 
	}

	int get_negsel_genotype(vector<int>& sites) { // returns number of seqs with a deleterious allele
		int geno = 0;		
		for (int i=0; i<2; ++i) {
			vector<int> intersectionality;
			set_intersection(sites.begin(), sites.end(), sequences[i].begin(), sequences[i].end(), back_inserter(intersectionality)); 
			if (intersectionality.size() >0) 
				++geno;
		}
		return geno;
	}
	// end new ch 8.4.1

	Individual (vector<vector<int>> seqs): sequences(seqs) {  // generation 0 and migration constructor
			;
	}

	Individual (Individual *p1, Individual *p2, vector<vector<int> > mutation_results, vector<int> breakpoints) {  // intra-simulation constructor //// NEW  CH6
		sequences.push_back((*p1).get_sequence(mutation_results[2][0]));    
		sequences.push_back((*p2).get_sequence(mutation_results[3][0]));    
		for (int i=0; i<2; ++i) {
			for (int j=1; j<mutation_results[i].size(); ++j)  {
				if (mutation_results[i][j] > 0) {
					sequences[i].push_back(mutation_results[i][j]);
					sort(sequences[i].begin(), sequences[i].end());   ///NEW CHAPTER 6 !!!!!!!!!!!!!!!!!!						
				}
				else 
					remove_allele_by_position(i, -1 * mutation_results[i][j]);
			}
		}

		//////////////NEW CHAPTER 6 //////////
		if (breakpoints.size() > 0)
			resolve_crossover( breakpoints );
		//////////// END NEW CHAPTER 6 ////////////////
	}

	~Individual() {}
};

#endif 
