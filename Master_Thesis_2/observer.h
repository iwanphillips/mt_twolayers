#ifndef observer_h
#define observer_h

double sync_error( const matrixd &x )
{
    size_t n = x.size1();
    double delta = 0.0;
    
    for( size_t i=0 ; i<n ; ++i )
    {
        delta += (x(i,1) - x(i,0))*(x(i,1) - x(i,0));
    }
    
    delta = sqrt( delta ); // /= n;

    return delta ;
}

struct observer
{
    ostream&        m_outfile;
    const double    m_sim_time;
    string          m_output;
    
    double m_delta;
    size_t m_count;
    double m_lamb2;
    
    observer(ostream &out, const float &sim_time, double lamb2 = 0.0) :
        m_delta( 0.0 ) , m_count( 0 ) ,
        m_outfile( out ) , m_sim_time( sim_time ) , m_lamb2( lamb2 ){}
    
    void set_lamb2( double lamb2 ) { m_lamb2 = lamb2; } //??
    
    void operator()(const matrixd &x, double t)
    {
        if(t > 100){                    // m_sim_time-5
            //if(fmod(t, 0.1) < 0.01){
            ++m_count;
            double mean = sync_error( x );
            m_delta += mean; // *m_dt
                
            if(t > m_sim_time-0.7){
                for(size_t i=0; i < x.size1(); ++i){
                    m_outfile << m_lamb2 << '\t' << t << '\t' << i << '\t' << x(i,0) << '\t' << x(i,1) <<  '\t' << (x(i,1) - x(i,0))*(x(i,1) - x(i,0)) << '\n';
                }
            }
            //}
        }
    }
    
    double get_delta( void ) const { return ( m_count != 0 ) ? m_delta / double( m_count ) : 0.0 ; }
    
    void reset( void ) { m_delta = 0.0; m_count = 0; }
};

/*
Attempt 2
m_delta.resize(m_delta.size()+1);
m_delta.insert_element(m_delta.size()-1, mean);
double sum = accumulate(delta.begin(), delta.end(), 0.0);
double mean = sum / delta.size();
cout << "final result is = " << mean << '\n';
 */

#endif /* observer_h */
