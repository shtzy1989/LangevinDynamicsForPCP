#ifndef MD_MATH_MATH_H
#define MD_MATH_MATH_H

#include "cppheader.h"
#include "MOLCONST.h"

#ifdef SINGLEPRECISION
typedef float Real;
#define Real_MAX FLT_MAX
#else
typedef double Real;
#define Real_MAX DBL_MAX
#endif

class Math{
public:
    // return the module rather than remainder by simple %
    static int Modulo(int x, int N){ return (x % N + N) % N; };
    static size_t Modulo(size_t x, size_t N){ return (x % N + N) % N; };
    static double Modulo(double x, double N){ return fmod(fmod(x, N) + N, N); }; // could be faster using if
    static float Modulo(float x, float N){ return fmod(fmod(x, N) + N, N); }; // could be faster using if
    // return the squre root
    static Real Sqrt(Real x){ return std::sqrt(x); };
    // return the inversed squre root
    static Real InversedSqrt(Real x){ return 1.0 / std::sqrt(x); };
    // return the cos
    static Real Cos(Real x){ return std::cos(x); };
    // return the sin
    static Real Sin(Real x){ return std::sin(x); };
    // return the acos
    static Real Acos(Real x){ return std::acos(x); };
    // return the exponential
    static Real Exp(Real x){ return std::exp(x); };
    // return the abs
    static Real Abs(Real x){ return std::abs(x); };

    template <typename T> static int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    static void GaussianRandom(Real *out1, Real *out2){
        Real s, u, v;
        do{
            u = Real(rand() % RAND_MAX) / Real(RAND_MAX) * 2.0 - 1.0;
            v = Real(rand() % RAND_MAX) / Real(RAND_MAX) * 2.0 - 1.0;
            s = u * u + v * v;
        }while( s == 0 || s >= 1 );

        if( out1 ) *out1 = 2.0 * u * sqrt(-0.5 * log(s) / s);
        if( out2 ) *out2 = 2.0 * v * sqrt(-0.5 * log(s) / s);
    };

    static Real Erf(Real x){ return std::erf(x); };
    static Real Erfc(Real x){ return std::erfc(x); };
    static Real Gaussian(Real x, Real u, Real s){ return 1.0 / (sqrt(2.0 * PI) * s) * exp(-0.5 * (x - u) * (x - u) / (s * s)); };
    static Real GaussianDerivative(Real x, Real u, Real s){ return -1.0 / (sqrt(2.0 * PI) * s * s * s) * (x - u) * exp(-0.5 * (x - u) * (x - u) / (s * s)); };
    static Real Integrate(Real (*f)(Real x), Real lowBound, Real highBound, Real step){
        Real sum = 0.0;
        for(Real x=lowBound;x<=highBound;x+=step){
            Real x0 = f(x);
            Real x1 = f(x + step);
            sum += (x0 + x1) * step * 0.5;
        }
        return sum;
    }

    static double IntegrateGaussianPositive(double x, double sigma){
        const double s2 = sqrt(2.0);
        return Erf(x / (s2 * sigma));
    };
    static double ErfInv(double a){
        // https://stackoverflow.com/questions/27229371/inverse-error-function-in-c
        double p, r, t;
        t = fma(a, 0.0 - a, 1.0);
        t = log(t);

        if( fabs(t) > 6.125 ){
            p =             3.03697567e-10; //  0x1.4deb44p-32 
            p = fma (p, t,  2.93243101e-8); //  0x1.f7c9aep-26 
            p = fma (p, t,  1.22150334e-6); //  0x1.47e512p-20 
            p = fma (p, t,  2.84108955e-5); //  0x1.dca7dep-16 
            p = fma (p, t,  3.93552968e-4); //  0x1.9cab92p-12 
            p = fma (p, t,  3.02698812e-3); //  0x1.8cc0dep-9 
            p = fma (p, t,  4.83185798e-3); //  0x1.3ca920p-8 
            p = fma (p, t, -2.64646143e-1); // -0x1.0eff66p-2 
            p = fma (p, t,  8.40016484e-1); //  0x1.ae16a4p-1 
        } else { // maximum ulp error = 2.35002
            p =             5.43877832e-9;  //  0x1.75c000p-28 
            p = fma (p, t,  1.43285448e-7); //  0x1.33b402p-23 
            p = fma (p, t,  1.22774793e-6); //  0x1.499232p-20 
            p = fma (p, t,  1.12963626e-7); //  0x1.e52cd2p-24 
            p = fma (p, t, -5.61530760e-5); // -0x1.d70bd0p-15 
            p = fma (p, t, -1.47697632e-4); // -0x1.35be90p-13 
            p = fma (p, t,  2.31468678e-3); //  0x1.2f6400p-9 
            p = fma (p, t,  1.15392581e-2); //  0x1.7a1e50p-7 
            p = fma (p, t, -2.32015476e-1); // -0x1.db2aeep-3 
            p = fma (p, t,  8.86226892e-1); //  0x1.c5bf88p-1 
        }
        r = a * p;
        return r;
    }

    static int32_t DivideRoundUp(int32_t dividend, int32_t divisor) {
        return (dividend + divisor - 1) / divisor;
    }
    static int64_t DivideRoundUp(int64_t dividend, int64_t divisor) {
        return (dividend + divisor - 1) / divisor;
    }


};


#endif
