#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <random>
#include <boost/random.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

typedef boost::numeric::ublas::vector<double> vectord;
typedef boost::numeric::ublas::matrix<double> matrixd;
using namespace boost::numeric::odeint;
using namespace std;

#include "system.h"
#include "observer.h"

int main(int argc, char **argv)
{
    string file = string(argv[1]);

    boost::property_tree::ptree desc;
    boost::property_tree::json_parser::read_json(file, desc);
    
    const int   n         = desc.get<int>("number of oscillators");
    const int   r_one     = desc.get<int>("range");
    const float o_one     = desc.get<float>("omega");
    const float g         = desc.get<float>("g");
    const float l_one     = desc.get<float>("lambda one");
    const float l_two     = desc.get<float>("lambda two");
    const float a     = desc.get<float>("alpha");
    const double sim_time = desc.get<float>("simulation time");
    string      path      = desc.get<string>("path");
    
    cout << "Num. Osc: "     << n      << '\n' << "Range: "               << r_one    << '\n' <<
            "Lambda one: "   << l_one  << '\n' << "Interlayer strength: " << l_two    << '\n' <<
            "Omega center: " << o_one  << '\n' << "Omega width: "         << g        << '\n' <<
            "Alpha: "        << a      << '\n' << "Simulation time: "     << sim_time << '\n';
     
    clock_t tStart = clock();
    uniform_real_distribution<double> distribution(-0.5,0.5);
    random_device rd;
    default_random_engine generator(rd());
    
    /*
    const size_t n = 256; int r_one = 35; float o_one = 0; float g = 0;
    float l_one = 0.085; float l_two = 0.045; float a = 1.45; const double sim_time = 100.0;
    string path = "/Users/iwanphillips/C++/Master_Thesis_Testing/2layer_data2.txt";
    */
    
    double dt = 0.025;
    matrixd x( n , 2 , 0. );
    
    k_ring system(n, l_one, l_two, o_one, g, r_one, a);
    
    string output_data = path;
    ofstream data_out(output_data);
    
    observer obs(data_out, sim_time, l_two);
    
    cout << "Starting integration..." << endl;
    
    //for( double lamb2 = 0.0 ; lamb2 < 0.1 ; lamb2 += 0.01){
    //for( double lamb1 = 0.0 ; lamb1 < 0.1 ; lamb2 += 0.01){
    
    for( double lamb2 = 0.0 ; lamb2 < 0.14 ; lamb2 += 0.005)
    {
        system.set_lamb2( lamb2 );
        obs.set_lamb2( lamb2 );
        obs.reset();
        
        for(size_t i = 0 ; i < n ; ++i )
        {
            double pos = i*2*M_PI/(n-1) - M_PI;
            double r1 = distribution(generator); double r2 = distribution(generator);
            x(i,0) = 6*r1*exp(-0.76*pos*pos); x(i,1) = 6*r2*exp(-0.76*pos*pos);
        }
        
        adams_bashforth_moulton<5, matrixd> abm;
        integrate_const(abm, system, x, 0.0, sim_time, dt, boost::ref( obs ));
        
        cout << lamb2 << "\t" << obs.get_delta() << endl;
    }

    data_out.close();
    
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
    
}

