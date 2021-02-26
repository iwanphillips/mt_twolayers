#ifndef observer_h
#define observer_h

typedef boost::numeric::ublas::matrix_column< matrixd > matrix_column;
#include <boost/numeric/ublas/matrix_proxy.hpp>

pair< double , double > calc_mean_field( const vectord &x , int range = 0 , int i = 1 )
{
    size_t n = x.size();
    double cos_sum = 0.0 , sin_sum = 0.0;
     
    for( int j=0 ; j<n ; ++j){
        float dist = abs(j-i);
        dist = abs(dist - round(dist/( (float) n ) ) * n);
             
        if(dist <= range && dist > 0){
        cos_sum += cos( x[j] );
        sin_sum += sin( x[j] );
        }
    }
    cos_sum /= double( 2 * (float) range );
    sin_sum /= double( 2 * (float) range );
         
     
    double K = sqrt( cos_sum * cos_sum + sin_sum * sin_sum );
    double Theta = atan2( abs(sin_sum) , abs(cos_sum) );
  
    return make_pair( K , Theta );
}

/* double sync_error( const matrixd &x )
{
    size_t n = x.size1();
    double psi_sum = 0.0;
    
    for( size_t i=0 ; i<n ; ++i )
    {
        double psi = fmod(abs(x(i,1)-x(i,0)), 2*M_PI);
        if(psi > M_PI){psi = 2*M_PI - psi;}
        psi_sum += psi * psi;
    }
    
    psi_sum /= n; // sqrt( psi_sum );

    return psi_sum;
} */

double om ( const matrixd &dxdt )
{
    size_t n = dxdt.size1();
    double psi_sum = 0.0;
    
    for( size_t i=0 ; i<n ; ++i )  // I may need modulus or similar here
    {
        double psi = dxdt(i,1)-dxdt(i,0);
        psi_sum += psi * psi;
    }
    
    psi_sum /= n;

    return psi_sum;
}

template<class system> // not necassary
struct observer
{
    ostream&        m_outfile;
    const double    m_sim_time;
    string          m_output;
    
    double m_om; // m_delta
    size_t m_count;
    double m_lamb2;
    double m_D;
    
    matrixd m_dxdt;
    system m_odefun;
    
    observer( system odefun, matrixd x, ostream &out, const float &sim_time, double lamb2 = 0.0 ) :
        m_odefun(odefun) , m_dxdt(x) , m_om( 0.0 ) , m_count( 0 ) ,
        m_outfile( out ) , m_sim_time( sim_time ) , m_lamb2( lamb2 ) { }
    
    void set_params( double lamb2 , double D ) { m_lamb2 = lamb2; m_D = D; }
    
    void operator()( matrixd &x, double t)
    {
        m_odefun( x, m_dxdt, t );
        
        cout << "t  = " << t << '\n';
        
        if(t > m_sim_time-5){ // 50
            ++m_count;
            
            //double mean = sync_error( x );
            double der_mean = om ( m_dxdt );
            
            m_om += der_mean;
                
            if(t > m_sim_time-0.05){
                for(int i=0; i < x.size1(); ++i){
                    
                    //matrix_column x1 (x,0); matrix_column x2 (x,1);
                    //pair< double , double > mean1 = calc_mean_field( x1 , 35 , i );
                    //pair< double , double > mean2 = calc_mean_field( x2 , 35 , i );
                    
                    m_outfile << m_D << '\t' << m_lamb2 << '\t' << t << '\t' << i << '\t' << // x(i,0) << '\t' << x(i,1) // << '\t' << (x(i,1) - x(i,0))*(x(i,1) - x(i,0))
                    // << '\t' << mean1.first << '\t' << mean1.second << '\t' << mean2.first << '\t' << mean2.second
                    m_dxdt(i,0) << '\t' << m_dxdt(i,1) << '\t' << m_dxdt(i,1) - m_dxdt(i,0) << '\n';
                    //cout << m_D << '\t' << m_lamb2 << '\t' << t << '\t' << i << '\t' << x(i,0) << '\t' << x(i,1) << '\t' << m_dxdt(i,1) << '\t' << m_dxdt(i,0) << '\n';
                    
                    //m_outfile << m_lamb2 << '\t' << t << '\t' << i << '\t' << m_dxdt(i,1) - m_dxdt(i,0) << '\n';
                }
            }
        }
    }
    
    double get_om( void ) const { return ( m_count != 0 ) ? sqrt( m_om / double( m_count ) ) : 0.0 ; } // sqrt in right place?
    
    void reset( void ) { m_om = 0.0; m_count = 0; }
};


#endif /* observer_h */
