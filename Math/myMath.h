//
//

#ifndef SLICEDOPTIM_MYMATH_H
#define SLICEDOPTIM_MYMATH_H

#include <functional>

long int binomialCoef(int n, int k);
long int binomialCoefNext(int n, int k, long int prev);
double myfmod(double v, double m);
double nballVolume(double r, int dim);
double nsphereArea(double r, int dim);

double integrateCosPower(double a, double b, int power);
double integrateSinPower(double a, double b, int power);

double simpsonIntegration(std::function<double(double)> f, double a, double b, int n);

#ifndef M_PI
#define M_PI 3.14159265358979323856
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#ifndef M_PI_4
#define M_PI_4 0.785398163397448309616
#endif
#ifndef M_1_PI
#define M_1_PI 0.318309886183790671538
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif


#endif //SLICEDOPTIM_MYMATH_H
