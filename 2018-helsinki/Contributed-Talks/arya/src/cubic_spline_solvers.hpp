var build_output_getCubicHittingTime(const double& t,
                 const std::vector<var>& theta,
                 const std::vector<double>& dtdtheta) {
  return precomputed_gradients(t, theta, dtdtheta);
}

double build_output_getCubicHittingTime(const double& t,
                    const std::vector<double>& theta,
                    const std::vector<double>& dtdtheta) {
  return t;
}

template <typename T0__, typename T1__, typename T2__, typename T3__, typename T5__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T5__>::type>::type
  getCubicHittingTime(const std::vector<T0__>& theta,
                      const T1__& t0,
                      const T2__& t1,
                      const T3__& threshold,
                      const int& NMAX,
                      const T5__& TOL, std::ostream* pstream__) {
    
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T5__>::type>::type RetType;
    // extract spline coefficients
    
    std::vector<double> theta_dbl = value_of(theta);
    double M0 = theta_dbl[0];
    double M1 = theta_dbl[1];
    double C0 = theta_dbl[2];
    double Ct0 = theta_dbl[3];
    
    std::vector<double> dtdtheta(theta.size());
    
    // bisection method
    int n = 1;
    double h = t1 - t0;
    double a = t0;
    double b = t1;
    double t, s_t, s_a;
    double dsdt;
    
    while(n <= NMAX) {
      
      t = (a+b)/2;
      s_t = (M0/(6*h))*pow(t1 - t, 3) + (M1/(6*h))*pow(t - t0, 3) + C0*(t-t0) + Ct0 - threshold;
      
      if(fabs(s_t) < TOL) {
        
        dsdt = -(M0/(2*h))*pow(t1 - t, 2) + (M1/(2*h))*pow(t - t0, 2) + C0;        
        dtdtheta[0] = -(1/(6*h))*pow(t1 - t, 3) / dsdt;
        dtdtheta[1] = -(1/(6*h))*pow(t - t0, 3) / dsdt;
        dtdtheta[2] = -(t-t0) / dsdt;
        dtdtheta[3] = -1.0 / dsdt;

        return  build_output_getCubicHittingTime(t, theta, dtdtheta);
      }
      
      s_a = (M0/(6*h))*pow(t1 - a, 3) + (M1/(6*h))*pow(a - t0, 3) + C0*(a-t0) + Ct0 - threshold;
      if((s_a > 0 && s_t > 0) || (s_a < 0 && s_t < 0)) {
        // a and t have the same sign
        a = t;
      } else {
        b = t;
      }
      n = n + 1;
    }
    
    // if max num iterations reached return error code
    return RetType(-42.0);
  }

var build_output_getMax(const double& m,
                                     const std::vector<var>& theta,
                                     const std::vector<double>& dmdtheta) {
  return precomputed_gradients(m, theta, dmdtheta);
}

double build_output_getMax(const double& m,
                                        const std::vector<double>& theta,
                                        const std::vector<double>& dmdtheta) {
  return m;
}

template <typename T0__, typename T1__, typename T2__, typename T4__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T4__>::type
  getMax(const std::vector<T0__>& theta,
         const T1__& t0,
         const T2__& t1,
         const int& NMAX,
         const T4__& TOL, std::ostream* pstream__) {
    
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T4__>::type RetType;
    
    // extract spline coefficients
    std::vector<double> theta_dbl = value_of(theta);
    double M0 = theta_dbl[0];
    double M1 = theta_dbl[1];
    double C0 = theta_dbl[2];
    double Ct0 = theta_dbl[3];
    
    std::vector<double> dmdtheta(theta.size());
    
    // bisection method
    int n = 1;
    double h = t1 - t0;
    double a = t0;
    double b = t1;
    double t, s_t, s_a;
    double max_s;
    
    while(n <= NMAX) {
      
      t = (a+b)/2;
      s_t = -(M0/(2*h))*pow(t1 - t, 2) + (M1/(2*h))*pow(t - t0, 2) + C0;
      //std::cout << "t: " << t << " s(t): " << s_t << std::endl;
      
      if(fabs(s_t) < TOL) {
        
        dmdtheta[0] = (1/(6*h))*pow(t1 - t, 3);
        dmdtheta[1] = (1/(6*h))*pow(t - t0, 3);
        dmdtheta[2] = (t-t0);
        dmdtheta[3] = 1.0;
        max_s = (M0/(6*h))*pow(t1 - t, 3) + (M1/(6*h))*pow(t - t0, 3) + C0*(t-t0) + Ct0;
        
        return  build_output_getMax(max_s, theta, dmdtheta);
      }
      
      s_a = -(M0/(2*h))*pow(t1 - a, 2) + (M1/(2*h))*pow(a - t0, 2) + C0;
      if((s_a > 0 && s_t > 0) || (s_a < 0 && s_t < 0)) {
        // a and t have the same sign
        a = t;
      } else {
        b = t;
      }
      n = n + 1;
    }
    
    // if max num iterations reached return error code
    return RetType(-42.0);
  }