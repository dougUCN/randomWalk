#ifndef OUTPUT
#define OUTPUT

#include <map>
#include "hdf5.h"
#include "hdf5_hl.h"
#include <iostream>
#include <vector>
#include <cstring>
#include <string>

const hsize_t NFIELDS = 10; // Change this when adding or removing columns from hdfOutputFormat
const int COMPRESSION = 0;
const int CHAR_COUNT = 16; // number of characters allowed for status field

struct hdfOutputFormat
{
     int    particleNum;
     int    windowHits;
     int    totalSteps;
     double location;
     char   status[CHAR_COUNT];
     int    cellRejections;
     int    cellExits;
     double velocity;
     double tstart;
     double tend;
     double mfp;

     hdfOutputFormat (  int pn, int wh, int ts, double loc,
                        std::string stringStatus, int cr,
                        int ce, double v, double tstart, double tend)
                        : particleNum(pn), windowHits(wh), totalSteps(ts), location(loc),
                        cellRejections(cr), cellExits(ce), velocity(v), tstart(tstart), tend(tend)
        {
            strncpy(status, stringStatus.c_str(), CHAR_COUNT); //need char[] for hdf c library to work
        }

};

void writeToHDF( std::string filename, std::vector<hdfOutputFormat> results , std::map<std::string, double> attributes);

// Print map to standard output
template<typename K, typename V>
void print_map(std::map<K,V> const &m)
{
    for (auto it = m.cbegin(); it != m.cend(); ++it) {
        std::cout << "{" << (*it).first << ": " << (*it).second << "}\n";
    }
}



#endif
