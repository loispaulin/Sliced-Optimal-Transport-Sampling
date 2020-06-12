//
//

#include <cmath>
#include <iostream>
#include <gsl/gsl_sf_hyperg.h>
#include "../Tools/my_utility.h"
#include "../Math/myMath.h"
#include "NBallRadonManager.h"

using namespace std;


/**
 * Compute the normalisation factor for the CDF of a \p d dimensional ball Radon transform
 * This amount to the ratio of the volume of a \p d - 1 dimensional unit ball over the one of a \p d dimensional unit ball
 * See wikipedia: Volume of an n-ball
 * @param d
 * @return
 */
double NBallRadonManager::computeNormalisationFactor(int d){
    if (d % 2){
        double res = d / 2.;
        for (int i = 1; i < (d+1)/2; ++i){
            res *= ((i - 0.5) / i);
        }
        return res;

    } else {
        double res = 1;
        for (int i = 1; i <= d/2; ++i){
            res *= (i / (i - 0.5));
        }
        res *= M_1_PI;
        return res;
    }
}

double NBallRadonManager::radon(double x) const{
    if (x < -1 || x > 1)
        return 0;
    double res = 0;
    double act = 1;
    for (int i = 0; i < (N - 1) / 2 + 1; ++i){
        res += act * coefs[i];
        act *= x*x;
    }
    if (N % 2){
        return res * normFactor;
    } else {
        return (res * sqrt(1 - x*x)) * normFactor;
    }
}

double NBallRadonManager::radonCDF(double x) const{
    if (x < -1)
        return 0.;
    if (x > 1)
        return 1.;
    if (N % 2){
        double res = 0;
        double act = x;
        for (int i = 0; i < (N - 1) / 2 + 1; ++i){
            res += act * CDFcoefs[i];
            act *= x*x;
        }
        return  res * normFactor - CDFoffset;
    } else {
        double res = std::sqrt(M_PI) * tgamma((1. + N) / 2.) / (2 * tgamma(1. + N / 2.));
        res += x * gsl_sf_hyperg_2F1(0.5, (1. - N) / 2., 1.5, x*x);
        return res * normFactor - CDFoffset;
    }
}


double NBallRadonManager::inverseCDF(double v) const{
    if (v >= 1)
        return 1;
    if (v <= 0)
        return -1;

    const double eps = std::pow(10, -10);

    double x = 0.;
    double w = radonCDF(x);

    while (abs(v - w) > eps) {
        double g = radon(x);
        double s = (v - w) / g;
        x += s;
        x = clamp(x, -1 + eps, 1. - eps);
        if (abs(s) < eps){
            break;
        }
        w = radonCDF(x);
        if (abs(x - 1.) < eps) {
            std::cerr << "x too close to one during cdf inversion" << std::endl;
            exit(13);
        }
    }

    return x;
}


NBallRadonManager::NBallRadonManager(int dimension){
    N = dimension;
    CDFoffset = 0;
    normFactor = computeNormalisationFactor(N);
    coefs.resize(N / 2 + 1);
    CDFcoefs.resize(N / 2 + 1);
    for (int i = 0; i < (N - 1) / 2 + 1; ++i){
        coefs[i] = binomialCoef((N - 1) / 2, i);
        if (i % 2){
            coefs[i] *= -1;
        }
        CDFcoefs[i] = coefs[i] / (2 * i + 1);

    }
    CDFoffset = radonCDF(-1);
}

