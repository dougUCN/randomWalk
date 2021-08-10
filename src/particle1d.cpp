#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <iostream>
#include <map>
#include "particle1d.hpp"

particle1d::particle1d( double startTime, double startVelocity, std::map<std::string, double> p )
            : start( p["start"] ), window( p["window"] ), cell( p["cell"] ),
            source( p["source"] ), cellChance( p["cellChance"] ), lossPerStep( p["lossPerStep"] ),
            windowLoss( p["windowLoss"] ), stepSize( p["stepSize"] ), stepSize2( p["stepSize2"] ) , fillTime( p["fillTime"] )
{
    resetState( startTime, startVelocity );
    sourceLeftCellRight = (source < cell);
}

void particle1d::resetState(  double startTime, double startVelocity, int num )
{
    particleNum = num;
    windowHits = 0;
    totalSteps = 0;
    cellExits = 0;
    cellEntranceChance = cellChance;
    location = start;
    prevLocation = start;
    status = "alive";
    cellRejections = 0;
    tstart = startTime;
    t = startTime;
    v = startVelocity;
    mfp = stepSize;
    cellExitLifetime = 4 / v * CELL_VOLUME / (std::pow(CELL_ENTRANCE_ID/2, 2 ) * M_PI);
}

std::map<std::string, boost::variant<double, int, std::string> > particle1d::getState()
{
    std::map<std::string, boost::variant<double, int, std::string>> output = {
        {"particleNum", particleNum},
        {"windowHits", windowHits},
        {"totalSteps", totalSteps},
        {"location", location},
        {"status", status},
        {"cellRejections", cellRejections},
        {"cellExits", cellExits},
        {"velocity", v},
        {"tstart", tstart},
        {"tend", t},
        {"mfp", t}
    };
    return output;
}

float particle1d::getLocation()
{
    return location;
}

void particle1d::walk( TMCGenerator &mc )
{
    while ( (t < fillTime) && status=="alive" )
    {
        // std::cout << totalSteps << "\tt: " << t << "\tloc: " << location << "\n";
        step( mc );
    }

}

void particle1d::step( TMCGenerator &mc )
{
    prevLocation = location;
    location +=  ( 1-(2*zeroOrOne(mc)) )  * mfp;
    t += mfp / v;
    totalSteps++;

    // Check collisions
    if (crossedWindow()) {
        windowHits++;
        if (zeroToOne(mc) < windowLoss) status = "window";  //chance for loss on the window
    } else if ( (location <= source && sourceLeftCellRight) || (location >= source && !sourceLeftCellRight) ) {
        status = "source"; // neutrons get absorbed by the source
    } else if ( (location >= cell && sourceLeftCellRight) ||  (location <= cell && !sourceLeftCellRight) ) {
        // chance for neutrons to get into cell
        if (zeroToOne(mc) < cellEntranceChance) {
            if ( (fillTime - t) < cellExitLifetime )
            {
                // TODO: option for cellExitLifetime to be calculated from exponential curve
                status = "cell";
                t = fillTime;
            } else {
                // If neutron exits the cell
                t += cellExitLifetime;
                cellExits++;
                totalSteps++;

                // Change MFP after cell exit
                if (stepSize2 != 0)
                {
                    mfp = stepSize2;
                    cellEntranceChance = std::pow( (1 - std::pow( 1- std::pow(CELL_ENTRANCE_ID/2 , 2) / (3 * 0.0254) / stepSize2 , 20) ) , 2 );

                }

                location = cell + ( mfp * leavingCellDirection());
            }
        } else {
        //neutron rejected from entrance
            location = prevLocation;
            totalSteps++;
            cellRejections++;
            t += mfp / v;
        }
    } else if (zeroToOne(mc) < lossPerStep ) {
        //chance to be absorbed by pipe
        status = "pipe";
    }


}

int particle1d::leavingCellDirection()
{
    if (sourceLeftCellRight) {
        return -1;
    } else {
        return 1;
    }
}

bool particle1d::crossedWindow()
{
    if (totalSteps <= 1) {
        // Ignore the first step
        return false;
    } else if (prevLocation < window &&  window < location) {
        return true;
    } else if (prevLocation > window  && window > location) {
        return true;
    } else if (location == window) {
        return true;
    } else {
        return false;
    }
}
