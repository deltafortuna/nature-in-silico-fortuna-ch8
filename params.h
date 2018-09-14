#ifndef PARAMS_H
#define PARAMS_H

// GLOBAL PARAMETERS
extern const int bitlength = 4000;  // must be constant at compile time; 
									//if you need a different bitlength size, you must recompile
									// Consider: After 10 genes, S was 2,394 for Ne=10,000, samplesize=500, r=1e-08, mu=1e-08
									// In this case, crashed without any error when bitlength set to 1000

extern double mutrate;
extern double recrate;    //////////// NEW CHAPTER 6 ????/////////
extern double hotrecrate; ////////////// NEW chapter 6 ///////
extern bool useRec;   //////////////// NEW CHAPTER 6 //////
extern bool useHotRec; ///////////// NEW Chapter 6 //////
extern bool trackAlleleBirths; //// NEW Chapter 7.4 //// 
extern int hotrecStart; /// NEw ch6 ///
extern int hotrecStop;  //// New ch6 ///
extern int seqlength;
extern int sampsize;
extern int sampfreq;
extern int getWindowStats; /////// new chp6 ///////
extern int windowSize; ///////////////// new chp 6 ////
extern int windowStep;  ///////// new chp 6   ///////
extern int pop_num; ///////////// new chp 7 //////////////
extern int runlength; //// new CHP 7 //////
extern Matrix<double> mig; ////////////////new chp 7 ////////////////
extern int diploid_sample; //// new ch 7.3  ///////
extern int printhapfreq; /// new ch 7.3
extern bool modelMigration; /// early7  //// 


// DEME=SPECIFIC PARAMETERS
int popsize; 
vector<int> demography;  // how to make all of these EXCEPT pop_schedule non-extern?????
vector<double> dem_parameter;
vector<int> dem_start_gen;
vector<int> dem_end_gen;
vector<int> carrying_cap;

extern vector<int> birthgen;
extern vector<int> extinctgen;
extern map<int, vector<int> > pop_schedule;
extern vector<bool> useMS; 
extern vector<string> mscommand;
extern map<int, vector<int> > splitgenesis;   //new ch 7.3
extern map<int, vector<int> > mergegenesis;   /// new chp 7.3
extern map<int, vector<double> > sellocus, possel, nfdsel, negsel; // new ch 8.4.1

#endif
