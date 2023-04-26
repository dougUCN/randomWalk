/**
* 1D random walk with a time element
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <map>
#include <iterator>
#include <boost/program_options.hpp>
#include <particle1d.hpp>
#include <boost/variant.hpp>
#include <boost/progress.hpp>
#include <output.hpp>
#include <chrono>
#include <mc.hpp>
#include <random>
#include <string>

namespace po = boost::program_options;

po::variables_map processArguments(int argc, const char** argv);

// uint64_t seed = 0; ///< random seed used for random-number generator (generated from high-resolution clock)

////////// input parameters ////////////
const double pipeID = 3 * 0.0254;                                // Inner diameter of the pipe [m]
const double pipeL  = 12;                                        // Total length of the pipe [m]
const double gateValve  = 6;                                     // starting position [m]
const double window = 9.1;                                       // PPM window location
const double source = 0;                                         // source location
const double fillTime = 50;                                    // source active time
////////////////////////////////////////

int main(int argc, const char *argv[])
{
    // Process cmd line args
    po::variables_map vm;
    try
    {
        vm = processArguments(argc,  argv);
    } catch (std::exception& err) {
        std::cerr << err.what() << "\nRun with --help argument\n";
        return 1;
    } catch (char const* helpFlag ) {
        std::cerr << helpFlag << '\n';
        return -1;
    }

     // Calculate parameters to pass to particle1d
    const double nonspec = vm["ns"].as<double>();                     // Chance for a nonspecular bounce
    const double lossPerBounce = vm["lpb"].as<double>();              // Loss per bounce (NiPh)
    const double windowLoss = vm["wl"].as<double>();                  // Chance for loss when passing through window
    const double mfp = pipeID * sqrt( 2*(2-nonspec)/nonspec/3 );     // Mean free path. eq.  4.79, eq.  4.70, and eq.  4.48 in Golub
    const double mfp2 =  vm["mfp2"].as<double>();                     // Mean free path after neutron exits cell
    // const double mfp = vm["mfp"].as<double>();     // Mean free path

    // get high-resolution timestamp to generate seed
    uint64_t seed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    // simulation time counter
    std::chrono::time_point<std::chrono::steady_clock> simstart = std::chrono::steady_clock::now();

    // load random number generator
    TMCGenerator mc;
    mc.seed(seed);
    std::uniform_real_distribution<double> start_time_distribution(0, nextafter(fillTime, std::numeric_limits<double>::max()) ); // Uniform distribution [0,fillTime]

    std::map<std::string, double> static_parameters = {
        {"stepSize",  mfp},
        {"nonspec",  nonspec},
        {"lossPerBounce",  lossPerBounce},
        {"cell", pipeL},
        {"window", window},
        {"windowLoss", windowLoss},
        {"start", gateValve},
        {"source", source},
        {"fillTime", fillTime},
        // {"cellChance", 0.35},
        {"cellChance", std::pow( (1 - std::pow( 1- std::pow(CELL_ENTRANCE_ID/2 , 2) / pipeID / mfp , 1/nonspec) ) , 2 )},
        {"stepSize2", mfp2},
        {"lossPerStep", (1/nonspec) * lossPerBounce},
    };

    std::cout << "\n### Parameters ###\n";
    print_map(static_parameters);

    std::vector<hdfOutputFormat> data;

    // Initialize particle and run through MC simulation
    particle1d ucn( start_time_distribution(mc), v2_average, static_parameters );
    boost::progress_display show_progress( vm["n"].as<int>()  );
    for  (int i = 0; i < vm["n"].as<int>(); i++)
    {
        ucn.walk( mc );
        std::map<std::string, boost::variant<double, int, std::string>> endstate = ucn.getState();
        // print_map( endstate );
        ucn.resetState( start_time_distribution(mc), v2_average, i + 1 );

        data.push_back( hdfOutputFormat{    boost::get<int>(endstate["particleNum"]),
                                            boost::get<int>(endstate["windowHits"]),
                                            boost::get<int>(endstate["totalSteps"]),
                                            boost::get<double>(endstate["location"]),
                                            boost::get<std::string>(endstate["status"]),
                                            boost::get<int>(endstate["cellRejections"]),
                                            boost::get<int>(endstate["cellExits"]),
                                            boost::get<double>(endstate["velocity"]),
                                            boost::get<double>(endstate["tstart"]),
                                            boost::get<double>(endstate["tend"]) }
                        );

        if (vm["progress"].as<bool>()) ++show_progress;
    }

    std::chrono::time_point<std::chrono::steady_clock> simend = std::chrono::steady_clock::now();
    float SimulationTime = std::chrono::duration_cast<std::chrono::milliseconds>(simend - simstart).count()/1000.;
    printf("\nSimulation Time: %.2fs\n", SimulationTime);

    // Output to file
    std::cout << "Writing data to file " << vm["f"].as<std::string>() << " ...";
    writeToHDF( vm["f"].as<std::string>(),  data, static_parameters);
    std::cout << "Done!\n";

    return 0;
}

/* Parse command line arguments */
po::variables_map processArguments(int argc, const char** argv)
{

    po::variables_map vm;
    po::options_description desc{"Usage"};
    desc.add_options()
        ("help,h", "Help")
        ("n", po::value<int>()->required(), "Number of particles to simulate")
        ("f", po::value<std::string>()->required(), "Stores end particle states to h5 file")
        ("lpb", po::value<double>()->default_value(1E-4), "Loss per bounce")
        ("ns", po::value<double>()->default_value(0.05), "Chance for nonspecular bounce")
        ("wl", po::value<double>()->default_value(0.03), "Chance of neutron loss for single window pass")
        // ("mfp", po::value<double>()->default_value(1), "1D walk step size")
        ("mfp2", po::value<double>()->default_value(0), "Mean free path upon cell exit")
        ("progress", po::value<bool>()->default_value(true), "Whether or not crude progress bar updates");;

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
        std::cout << desc << '\n';
        throw "Doug Wong 2021";
    } else {
        po::notify(vm); // Must be after vm.count("help")
    }

    return vm;
}
