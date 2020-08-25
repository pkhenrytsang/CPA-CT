#include<complex>

static inline double sqr(std::complex<double> x){ return sqr(std::abs(x)); }

static inline double sqr(double x){ return x*x; }

static inline double quar(double x){ return x*x*x*x; }


