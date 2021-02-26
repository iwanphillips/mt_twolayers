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

#include "noise.h"
#include "system.h"
#include "observer.h"

int main(int argc, char **argv)
{
    /*
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
     */
     
    clock_t tStart = clock();
    uniform_real_distribution<double> distribution(-0.5,0.5);
    random_device rd;
    default_random_engine generator(rd());
    
    const size_t n = 256; int r_one = 35; float o_one = 1; float g = 0;
    float l_one = 0.02; float l_two = 0.01; float a = 1.3908; const double sim_time = 400.0; // l_two = 0.045; // a = 1.45 // l_one = 0.085
    string path = "/Users/iwanphillips/C++/mt4_odeint/D0_l1_02_l2_0_01_gamma0p1.txt";
    double D = 0.0; double gamma = 0.1;
    
    double dt = 0.01; // 0.025;
    matrixd x( n , 2 , 0. );
    //double alpha_ = 1;
    
    //
        int steps = (float) sim_time / dt;
        int n_pts = 16384;
    /*
        while(steps > n_pts){n_pts = 2*n_pts;}
        //float Q_d = 1;
        
        int c = 5;
        int max = c*n_pts;
        while(max < 2000000){max += steps; c += 1;}
        vector<float> X(max, 0); // X //float X[max];
        for( int i=0 ; i<c ; ++i ){ //for( int i=0 ; i<50 ; ++i ){ //i<n
            float X_i[n_pts];
            memset(X_i, 0.0, sizeof X_i);
            long utime; utime=(long)time(NULL); long seed=utime;
            f_alpha(n_pts, X_i, Q_d, alpha_, &seed);
            vector<double> Xi (X_i, X_i + sizeof(X_i) / sizeof(int) );
            X.insert(X.end(), Xi.begin(), Xi.end());
        }
    */
       int c = 10; int max = c*steps;
       while(max < 2000000){max += steps; c += 1;}
       vector<float> X;
       for( int i=0 ; i<c ; ++i ){
           float X_i[steps];
           memset(X_i, 0.0, sizeof X_i);
           long utime; utime=(long)time(NULL); long seed=utime;
           cor_exp(steps, X_i, 0.01, &seed, dt, gamma);
           vector<double> Xi (X_i, X_i + sizeof(X_i) / sizeof(int) );
           X.insert(X.end(), Xi.begin(), Xi.end());
       }
    //
    
    k_ring system( X, n, l_one, l_two, o_one, g, r_one, a, D ); // not a function call, even though it looks like one
    
    string output_data = path;
    ofstream data_out(output_data);
    
    observer< k_ring > obs(system, x, data_out, sim_time, l_two);
    
    cout << "Starting integration..." << endl;
    
    for(size_t i = 0 ; i < n ; ++i )
    {
        double pos = i*2*M_PI/(n-1) - M_PI;
        double r1 = distribution(generator); double r2 = distribution(generator);
        x(i,0) = 6*r1*exp(-0.76*pos*pos); x(i,1) = 6*r2*exp(-0.76*pos*pos);
    }
    
    /*
    adams_bashforth_moulton<5, matrixd> abm;
    integrate_const(abm, system, x, 0.0, sim_time, dt, boost::ref( obs )); //runge_kutta4< matrixd >()

    data_out << D << '\t' << l_two << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 <<  '\t' << obs.get_om() << '\n';
    //
    for( double l_two = 0 ; l_two < 0.0101; l_two += 0.001 ){
        system.set_params( l_two , 0 );
        obs.set_params( l_two , 0  );
        obs.reset();
    
        adams_bashforth_moulton<5, matrixd> abm;
        integrate_const(abm, system, x, 0.0, sim_time, dt, boost::ref( obs )); //runge_kutta4< matrixd >()
    
        //cout     << D << '\t' << l_two << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 <<  '\t' << obs.get_om() << '\n';
        //data_out << D << '\t' << l_two << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 <<  '\t' << obs.get_om() << '\n';
    }
    */
     for( double l_two = 0 ; l_two < 0.050; l_two += 0.001){
        for( double D_perc  = 0.0 ; D_perc < 0.1001 ; D_perc += 0.01 ){
            system.set_params( l_two , D_perc*l_two  );
            obs.set_params( l_two , D_perc*l_two  );
            obs.reset();
        
            adams_bashforth_moulton<5, matrixd> abm;
            integrate_const(abm, system, x, 0.0, sim_time, dt, boost::ref( obs )); //runge_kutta4< matrixd >()
        
            cout     << D_perc << '\t' << l_two << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 <<  '\t' << obs.get_om() << '\n';
            data_out << D_perc << '\t' << l_two << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 <<  '\t' << obs.get_om() << '\n';
        }
    }
    //

    data_out.close();
    
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
    
}

