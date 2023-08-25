#include <iostream>
#include <omp.h>
#include "ShallowWater.h"

// add -lboost_program_options to compile
#include <boost/program_options.hpp>
namespace po = boost::program_options;



int main(int argc, char* argv[]) {

    int check;

    // Gets the number of threads, as defined in OMP_NUM_THREADS env variable.
    // int max_num_threads = omp_get_max_threads();
    
    // if (max_num_threads > 40) {
    //         std::cout << "Error: Program is being run with " << max_num_threads << " threads." << std::endl;
    //         std::cout << "The program must be run with 40 or less threads." << std::abort;
    //         std::cout << "Please change the number of threads by altering the 'OMP_NUM_THREADS' env variable." << std::endl;
    //         return 0;   // Exiting with 0 as error was properly handled.
    // }
    // else {  // Displaying on terminal how many threads will be used.
    //     std::cout << "Running program with " << max_num_threads << " thread(s)." << std::endl;
    // } 

    // declare variables for command line input (and default values)
    po::options_description opts("Allowed options");
    opts.add_options()
        ("h", "produce help message")

        ("dt", po::value<double>()->required(), "time-step (double)")

        ("T", po::value<int>()->required(), "total integration time (int.)")

        ("Nx", po::value<int>()->default_value(100), "number of grid points in x (int.)")

        ("Ny", po::value<int>()->default_value(100), "number of grid points in y (int.)")

        ("ic", po::value<int>()->default_value(4), "Index of initial condition (1-4)")

        ("choice", po::value<int>(&check)->default_value(1), "For loop (1) or Blas (2)");

    po::variables_map vm;

     // Handles errors with the passed arguments and provides helpful error messages.
    try {

        po::store(po::parse_command_line(argc, argv, opts), vm);
        // If the user asks for help, or no arguments are supplied, provide
        // a help message and list the possible arguments.
        if (vm.count("h") || argc == 1) {
            std::cout << "Please provide a value for all the required arguments " << std::endl
                    << opts << std::endl;
            return 0;   // Exiting with 0 as error was properly handled.
        }

        po::notify(vm);

        // check if the user input the correct option for solution methods
        if (check != 1 && check != 2)
        {
            std::cout << "Invalid solution method input" << std::endl
                        << opts << std::endl;
            return 1;
        }
        
    }
    catch (std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        std::cout << opts << std::endl;
        return 0;
    }

    // Extracts the parameter values given to variables using the appropriate datatype.
    const double dt = vm["dt"].as<double>();
    const int T = vm["T"].as<int>();
    const int Nx = vm["Nx"].as<int>();
    const int Ny = vm["Ny"].as<int>();
    const int ic = vm["ic"].as<int>();
    const int choice = vm["choice"].as<int>();




    ShallowWater sha_water;


    sha_water.SetParameters(dt, T, Nx, Ny, ic, choice);

    sha_water.SetInitialConditions();

    sha_water.TimeIntegrate();

    sha_water.SaveToFile();
    


    return 0;

}




