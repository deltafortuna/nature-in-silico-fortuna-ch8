# nature-in-silico-fortuna-ch8
FORTUNA simulation environment for chapter 8
The code in this directory contains FORTUNA files updated for additional functionality as 
described in Chapter 8 of the book. 

Note that the parameters file has additional entries as well. As usual, you will need a 
compiled copy of MS in the same directory as the one where you compile FORTUNA if you
are using MS to generate the starting genetic variation of one or more demes. 

To compile:
g++ -std=c++11 fortuna.cc -o fortuna
