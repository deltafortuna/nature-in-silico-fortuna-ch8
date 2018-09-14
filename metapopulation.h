#ifndef METAPOPULATION_H   /// New class starting with CH 7
#define METAPOPULATION_H
/*
using std::set_difference;
using std::inserter;*/

class Metapopulation {

private:
vector<Population*> populations;
uniform_real_distribution<double> random01;

public:

	int reproduce_and_migrate(int gen) {  // ch 8.4.1 changed from void to int. 
		bool fixtest = true; // new ch 8.4.1
		// reproduction within all extant demes and check for birth/extinction of population 
		for (int i=0; i<pop_num; ++i) {
			if ((*populations[i]).get_extant()) {
				if (extinctgen[i] == gen) 
					(*populations[i]).set_extinct();  // free memory???
				else
					(*populations[i]).reproduce(gen);
//cout << (*populations[i]).get_individuals_vectorsize() << " in population #" << i << " after reproduction." << endl;								
			} else {
//// alllow set_extant to check for split or merger AND to initialize individuals in new population if so. 				
				if ( birthgen[i] == gen ) {
					// new ch 7.3
					vector<int> change = (*populations[i]).set_extant();
					map<int, int> alleles_to_add;					
					if (change[0] == 1) { //split
						int N = (double) (change[2])  / 100 * (*populations[ change[1] ]).get_current_popsize(gen);  // does this round to lowest number
						for (int k = 0; k<N; ++k) {						
							vector<vector<int>> v1 = (*populations[ change[1] ]).get_sequences(k);
							(*populations[i]).add_immigrant(  v1  );		
							for (int m=0; m<2; ++m) 
								for (auto iter=v1[m].begin(); iter!=v1[m].end(); ++iter)
									++alleles_to_add[*iter];			
						}
						(*populations[ change[1] ]).remove_emigrants(N); // so that the other deme populated by the split doesn't receive same individuals			
					}	else if (change[0] == 2) { //merger
						vector<int> N;
						N.push_back( (*populations[ change[1] ]).get_current_popsize(gen) );					
						N.push_back( (*populations[ change[2] ]).get_current_popsize(gen) );
						for (auto q:N) {						
							for (int k = 0; k<q; ++k) {
								vector<vector<int> > v1 = (*populations[ change[1] ]).get_sequences(k);
								(*populations[i]).add_immigrant( v1 );
								for (int m=0; m<2; ++m) 
									for (auto iter=v1[m].begin(); iter!=v1[m].end(); ++iter)
										++alleles_to_add[*iter];
							}
						}	
					}	 // else built from MS	
					if (change[0] > 0) // deme not built from MS
						for (auto iter=alleles_to_add.begin(); iter != alleles_to_add.end(); ++iter) // populate alleles in new deme
							(*populations[i]).insert_new_allele( (*populations[ change[1] ]).get_allele_info(iter->first) );
					// end // new ch 7.3	
					(*populations[i]).reproduce(gen);  				
				}
			} 
		}

		// migration among all demes
		for (int i=0; i<pop_num; ++i) {
			if ( (*populations[i]).get_extant() ) {
				for (int j=0; j < pop_num; ++j) {	
					if ( (*populations[j]).get_extant() ) {
						//double Nm =0;		
						double Nm = mig[i][j] * pop_schedule[j][gen]; 
						if (Nm < 1) {
							if(random01(f) < Nm)
								Nm = 1; 
							else 
								Nm = 0;
						}
						else 
							Nm = floor(Nm); 

						for (int k=0; k<Nm; ++k) {				
							vector<vector<int>> v1 = (*populations[i]).get_sequences(k);
							(*populations[j]).add_immigrant(  v1  );

							// check for new alleles introduced to population j	
							for (int m = 0; m < 2; ++m) {			
								vector<int> v2 = (*populations[j]).get_allele_positions();
								vector<int> diff;
								set_difference(v1[m].begin(), v1[m].end(), v2.begin(), v2.end(), inserter(diff, diff.begin()));
								for (auto q:diff) 
										(*populations[j]).insert_new_allele( (*populations[i]).get_allele_info(q) ) ;												
							}
						}
// if (i == 1 && j==2) cout << "Gen " << gen <<  " [1,2], pop 1 is " << 						
						(*populations[i]).remove_emigrants(Nm);
					} else continue;
				} 
			} else continue;	
		} 
					
		//// get samples when appropriate   /// early7
		if ( gen==0 || (gen+1) % sampfreq == 0)  {
			for (int i=0; i<pop_num; ++i) {				
				if ((*populations[i]).get_extant()) {					
					if (possel[(*populations[i]).get_popnum()][0] != 0 && sellocus[(*populations[i]).get_popnum()][4] == 1 )
						fixtest = (*populations[i]).sample(gen);	
					else 
						int dummy = (*populations[i]).sample(gen);
					if (! fixtest )
						break;
				}
			}
		}		
		return (fixtest);	// new ch 8.4.1
	}

	void close_output_files() {
		for (auto iter = populations.begin(); iter != populations.end(); ++iter)
			(*iter)->close_output_files();
	}

	Metapopulation()  {
		for (int i=0; i < pop_num; ++i) {
			if (birthgen[i] != 0 ) 
				populations.push_back( new Population(i, 0) );
			else 
				populations.push_back( new Population(i, 1) );
		}
		random01.param(uniform_real_distribution<double>::param_type(0.,1.));
	}
	~Metapopulation() {}

	static mt19937 f;
};

#endif 
