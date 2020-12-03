#ifndef system_h
#define system_h

struct k_ring
{
    int     m_N;
    double  m_lamb1;
    double  m_lamb2;
    vectord m_omega;
    int     m_range;
    float   m_alpha;
    matrixd m_loc_coup;
    matrixd m_int_coup;
    
    k_ring( int N = 5 , double lamb1 = 1.0 , double lamb2 = 1.0 , float om_loc = 0.5, double g = 0.0 , int range = 1 , // double omega = 2.0
           float alpha = M_PI-0.1 ) :
    m_N ( N ) , m_lamb1 ( lamb1 ) , m_lamb2 ( lamb2 ) , m_omega ( N , 0.0 ) , m_range ( range ) , // m_omega ( omega )
    m_alpha ( alpha ) , m_loc_coup(N , 2, 0.) , m_int_coup(N , 2, 0.)
    {
        create_frequencies( om_loc, g );
    }
    
    void create_frequencies( float om_loc, double g )
    {
        boost::mt19937 rng;
        boost::cauchy_distribution<> cauchy( om_loc , g );
        boost::variate_generator< boost::mt19937&, boost::cauchy_distribution<> > gen( rng , cauchy );
        generate( m_omega.begin() , m_omega.end() , gen );
    }
    
    void set_lamb2( double lamb2 ) { m_lamb2 = lamb2; }

    double get_lamb2( void ) const { return m_lamb2; }
    
    void operator() (const matrixd &x, matrixd &dxdt, const double )
    {
        for( int i=0 ; i<m_N ; ++i ){
            
            m_int_coup(i,0) = sin( x(i,0) - x(i,1) );
            //m_int_coup(i,1) = sin( x(i,1) - x(i,0) );
            
            m_loc_coup(i,0) = m_loc_coup(i,1) = 0.;
            //double a = i*2*M_PI/(m_N-1) - M_PI;
            
            for( int k=0; k<m_N ; ++k ){
                //double b = k*2*M_PI/(m_N-1) - M_PI;
                float dist = abs(k-i);
                dist = abs(dist - round(dist/( (float) m_N ) ) * m_N);
                if(dist <= m_range && dist > 0){
                    //m_loc_coup(i,0) += (1/(2*M_PI))*(1 + 0.995*cos(a-b))*sin( x(i,0) - x(k,0) + m_alpha );
                    //m_loc_coup(i,1) += (1/(2*M_PI))*(1 + 0.995*cos(a-b))*sin( x(i,1) - x(k,1) + m_alpha );
                    m_loc_coup(i,0) += sin( x(i,0) - x(k,0) + m_alpha );
                    m_loc_coup(i,1) += sin( x(i,1) - x(k,1) + m_alpha );
                }
            }
            dxdt(i,0) = m_omega(i) - ((m_lamb1)/(2*m_range + 1))*m_loc_coup(i,0) + (m_lamb2/2)*m_int_coup(i,0); // (2*M_PI*m_loc_coup(i,0)/(double) m_N) // ((m_lamb1)/(2*m_range + 1))
            dxdt(i,1) = m_omega(i) - ((m_lamb1)/(2*m_range + 1))*m_loc_coup(i,1) - (m_lamb2/2)*m_int_coup(i,0); // (1/(2*M_PI))* ?? // /(double m_N) ?? //((m_lamb1)/(2*m_range + 1))*m_loc_coup(i,1)
        }
    }
};

#endif /* system_h */
