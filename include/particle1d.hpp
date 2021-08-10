#ifndef PARTICLE1D
#define PARTICLE1D

#include <map>
#include <boost/variant.hpp>
#include <mc.hpp>

const double CELL_ENTRANCE_ID = 0.0745; // [meters]
const double CELL_VOLUME = 0.019635;    // [m^3]
const double v2_average = 4.81342; // Average velocity of a v2 dv distribution up to 214 neV [m/s]
const double v3_average = 5.13431; // Average velocity of a v2 dv distribution up to 214 neV [m/s]
const double v2_average_SS = 4.44; // Average velocity of a v2 dv distribution up to 190 neV [m/s]
const double v3_average_SS = 4.73; // Average velocity of a v3 dv distribution up to 190 neV [m/s]

class particle1d {
public:
    particle1d( double startTime, double startVelocity, std::map<std::string, double> p);
    void resetState( double startTime, double startVelocity, int num = 0 );
    std::map<std::string, boost::variant<double, int, std::string> > getState();
    void walk( TMCGenerator &mc );
    float getLocation();

private:
    int particleNum;
    double t, v, location, prevLocation, tstart, mfp;
    double cellExitLifetime;            // Lifetime for a neutron at velocity v to exit the precession cell
    double cellEntranceChance;          // Chance for the neutron to enter the cell 
    std::string status;
    bool sourceLeftCellRight;
    int cellRejections, windowHits, totalSteps, cellExits;

    const double start;                  // Neutron start location
    const double window;                 // location of window
    const double cell;                   // location of cell
    const double source;                 // location of source
    const double cellChance;             // Chance for particle to enter cell
    const double lossPerStep;            // Loss per random walk step
    const double windowLoss;             // Single pass chance of window loss
    const double stepSize;               // 1D walk step size
    const double stepSize2;               // 1D walk step size after exiting from cell
    const double fillTime;               // Source active time

    void step( TMCGenerator &mc );                        // Takes a 1D step
    int leavingCellDirection();                // Returns -1 if sourceleftCellRight is true, +1 if false
    bool crossedWindow();                // Returns true if neutron crossed window

};


#endif
