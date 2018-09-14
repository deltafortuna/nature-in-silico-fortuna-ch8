#include "params.h"  // access to declaratiosns of global parameter values 

map<int, map<string, string> > read_parameters_file(const string &parameters_fn) 
{	
	map<int, map<string, string> > params_by_block;
	map<string,string> params;
	ifstream paramfile(parameters_fn.c_str());
	string line;
	int block = -1;  // -1 for global parameters, non-neg blocks are deme-specific parameters
	regex query("DEME");
	while(getline(paramfile, line)) {	
		istringstream iss(line);
		string key, nextone, value;
		iss >> key;
		if (regex_search(key, query)) {  // check if entering next deme block		
			params_by_block[block] = params;
			params.clear();
			++block;			
		} else { 
			while (iss >> nextone)
				value += nextone + " ";
			params[key] = value;
		}
	}		
	params_by_block[block] = params; // assign values for final deme
	return  params_by_block;
}

map<int, map<string, string> > p = read_parameters_file("parameters"); // change string literal argument if you change the parameters file name to something else

vector<int> get_multi_int_param(const string &key, map<string, string> &parameters)    /// USE TEMPLATE SO ONLY ONE OF THESE FUNCTIONS IS NEEDED
{
	vector<int> vec;
	istringstream iss(parameters[key].c_str());
	string param;
	while(getline(iss, param, ' '))
		vec.push_back(atoi(param.c_str()));
	return vec;
}

vector<double> get_multi_double_param(const string &key, map<string, string> &parameters)
{
	vector<double> vec;
	istringstream iss(parameters[key].c_str());
	string param;
	while(getline(iss, param, ' '))
		vec.push_back(atof(param.c_str()));
	return vec;
}

vector<int> create_pop_schedule() 
{
	vector<int> ps;
	int i=0;
	int cursize = popsize;
	for (int step = 0; step < demography.size(); ++step) {
		for (; i<dem_start_gen[step]; ++i) 
			ps.push_back(cursize);
		for (; i <= dem_end_gen[step]; ++i) {
			switch(demography[step]) {
				case 0: ps.push_back(cursize);  // no size change
					    break;
				case 1: if (i == dem_start_gen[step])
							cursize += dem_parameter[step]; 
						ps.push_back(cursize); // instantaneous change 
						break;  
				case 2: cursize += dem_parameter[step];
						ps.push_back(cursize);  // linear change
						break;
				case 3: cursize *= exp(dem_parameter[step]);  // exponential change
						ps.push_back(cursize);
						break;
				case 4: cursize = (carrying_cap[step] * cursize) / 
									(cursize + (carrying_cap[step] - cursize)*exp(-1*dem_parameter[step]));  // logistic change
						cout << i << " " << cursize << endl;
						ps.push_back(cursize);
			}
		}
	}
	return ps;
} 

/////////// END NEW ch 7 ////////
double mutrate, recrate, hotrecrate;
int hotrecStart, hotrecStop, sampsize, seqlength, sampfreq, getWindowStats, windowSize, windowStep, pop_num, runlength, diploid_sample, printhapfreq; // new ch 7.3
bool useRec, useHotRec, modelMigration, trackAlleleBirths;   // ch 7.4;   /// early7 ////
vector<double> m;
vector<string> mscommand;
vector<bool> useMS; 
vector<int> birthgen, extinctgen; 
map<int, vector<int> > pop_schedule, splitgenesis, mergegenesis; // new ch 7.3;
map<int, vector<double> > sellocus, possel, nfdsel, negsel; /// new ch 8.4.1 


int process_parameters() {
	for (auto iter = p.begin(); iter!=p.end(); ++iter ) {
		map<string, string> parameters = iter->second; 	
		if (iter->first == -1) {  // block of global parameters
			mutrate = atof(parameters["mutrate"].c_str());
			recrate = atof(parameters["recrate"].c_str());    ////////////////////?NEW chapter 6 /////
			hotrecrate = atof(parameters["hotrecrate"].c_str()); // New chapter 6 ///////////
			useRec = atoi(parameters["useRec"].c_str());  //// NEW chapter 6 /////
			useHotRec = atoi(parameters["useHotRec"].c_str()); /////////////NEW chapter 6 ////////
			hotrecStart = atoi(parameters["hotrecStart"].c_str());  //// New chp 6 /////
			hotrecStop = atoi(parameters["hotrecStop"].c_str());  //// New chp 6 /////
			sampsize = atoi(parameters["sampsize"].c_str());
			seqlength = atof(parameters["seqlength"].c_str()); // covernsion using atof() enables use of e notation in parameters file
			sampfreq = atoi(parameters["sampfreq"].c_str());  
			getWindowStats = atoi(parameters["getWindowStats"].c_str());
			windowSize = atoi(parameters["windowSize"].c_str());
			windowStep = atoi(parameters["windowStep"].c_str());			
			pop_num = atoi(parameters["pop_num"].c_str());        /////// nnew chapter 7 //////	
			runlength = atoi(parameters["runlength"].c_str());
			m = get_multi_double_param("migration_rates", parameters);   //// new chp 7 ////////	
			diploid_sample = atoi(parameters["diploid_sample"].c_str()); // /new ch 7.3 ////// 
			printhapfreq = atoi(parameters["printhapfreq"].c_str()); // new ch 7.3
			trackAlleleBirths = atoi(parameters["trackAlleleBirths"].c_str()); // new ch 7.4 ///
			modelMigration = atoi(parameters["modelMigration"].c_str()); /// early7 // 
		} else {   // block of deme parameters
			popsize = atoi(parameters["popsize"].c_str());
			demography =  get_multi_int_param("demography", parameters); 
			dem_parameter = get_multi_double_param("dem_parameter", parameters); 
			dem_start_gen =  get_multi_int_param("dem_start_gen", parameters); 
			dem_end_gen = get_multi_int_param("dem_end_gen", parameters); 
			carrying_cap = get_multi_int_param("carrying_cap", parameters); 
			birthgen.push_back( atoi(parameters["birthgen"].c_str()) );
			extinctgen.push_back( atoi(parameters["extinctgen"].c_str()) );
			useMS.push_back( atoi(parameters["useMS"].c_str()) );
			mscommand.push_back( parameters["mscommand"] );
			pop_schedule[iter->first] = create_pop_schedule();
			splitgenesis[iter->first] = get_multi_int_param("splitgenesis", parameters);		
			mergegenesis[iter->first] = get_multi_int_param("mergegenesis", parameters);
			sellocus[iter->first] = get_multi_double_param("sellocus", parameters); // new ch 8.4.1
			possel[iter->first] = get_multi_double_param("possel", parameters);  // new ch 8.4.1
			nfdsel[iter->first] = get_multi_double_param("nfdsel", parameters); // new ch 8.4.1
			negsel[iter->first] = get_multi_double_param("negsel", parameters); // new ch 8.4.2
		}
	}
	return 1;
}

int good_parameters = process_parameters();
double* a = &m[0]; // matrix takes array argument
Matrix<double> mig(pop_num, pop_num, a);