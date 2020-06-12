//
//

#ifndef SLICEDOPTIM_MAPPING_H
#define SLICEDOPTIM_MAPPING_H

#include "../Math/VecX.h"
#include "../Math/myMath.h"


//We have numerical instability when using dimension > 40. Using long double everywhere fixes it
double my_gamma(int dim);
double tau(int dim);
double rho(int dim);
double tau(double lambda, int dim);
double dtau(double lambda, int dim);
double inverseTau(double lambda, int dim);
double rho(double lambda, int dim);
double drho(double lambda, int dim);
double inverseRho(double lambda, int dim);
double solveForGamma(int dim);

/*
template<typename VecType >
inline double gtopgsl(const VecType& x, double y) {
    int d = int(x.dim());
    double y2 = y * y;
    double x2 = x.norm2();
    double y2_x2 = y2 / x2;
    double sqrt1 = std::sqrt(1 + y2_x2);

    double leftSub = (std::sqrt(M_PI) * tgamma(0.5* d)) / (2. * tgamma((1 + d) * 0.5));
    double h2F1;
    if (d == 2){
        h2F1 = 1;
    } else {
        h2F1 = gsl_sf_hyperg_2F1(0.5, 1. - 0.5 * d, 1.5, y2 / (x2 + y2));
    }
    double rightSub = (y * h2F1) / (x.norm() * sqrt1);

    return sqrt1 * std::pow(d * (leftSub - rightSub), 1. / d);
}
 */
template<typename VecType>
inline double gtop(const VecType& x, double y) {
    double x1 = x.norm();
    double x2 = x1 * x1;
    static const double eps = std::pow(10, -16);
    if (x2 < eps){
        return 1.;
    } else {
        double factor = std::sqrt((y*y / x2) + 1.);
        int d = x.dim();
        return factor * std::pow(d * integrateSinPower(0, std::atan(x1 / y), d-1), 1./d);
    }
}

template<typename VecType>
inline double gbot(const VecType& x, double y) {
    return std::sqrt(1 + y * y / x.norm2());
}

template<typename VecType >
inline double htop(const VecType& x, double y) {
    return std::sqrt(x.norm2() + y * y);
}
/*
template<typename VecType >
inline double hbotgsl(const VecType& x, double y) {
    int d = int(x.dim());
    double y2 = y * y;
    double x2 = x.norm2();
    double y2_x2 = y2 / x2;

    double leftSub = std::sqrt(M_PI) / std::tgamma((1. + d) * 0.5);
    double rightSub = gsl_sf_hyperg_2F1_renorm(0.5, 0.5 * d, 1. + 0.5 * d, 1. / (1. + y2_x2))
                      * std::pow(1. + y2_x2, -0.5 * d);

    return 0.5 * std::sqrt(y2 + x2) * tgamma(0.5 * d) * (leftSub - rightSub);
}
 */
template<typename VecType >
inline double hbot(const VecType& x, double y) {
    double x1 = x.norm();
    double x2 = x1 * x1;
    double y2 = y * y;
    double factor = std::sqrt(x2 + y2);
    int d = x.dim();
    return factor * integrateCosPower(0, std::atan(y / x1), d-1);

}

inline VecX<2> polar2carthesian(const VecX<2> &p){
    return VecX<2>(p[0] * cos(p[1]), p[0] * sin(p[1]));
}
inline VecX<2> square2disk(const VecX<2>& p){
    if (p.dim())
        if (p.norm() < std::pow(10, -16)){
            return p;
        }
    if(std::abs(p[0]) > std::abs(p[1])){
        return polar2carthesian(VecX<2>(p[0], p[1] / p[0] * M_PI_4));
    } else {
        return polar2carthesian(VecX<2>(p[1], (2 - p[0] / p[1]) * M_PI_4));
    }
}

template<int DIM>
inline VecX<DIM> NBall2NCylinder(const VecX<DIM>& p) {
    if (p.norm() < std::pow(10, -16)) {
        return p;
    }
    if (p.dim() == 1) {
        return p;
    }
    VecX<DIM-1> x;
    double y = p[DIM - 1];
    for (int i = 0; i < DIM - 1; ++i) {
        x[i] = p[i];
    }

    double signy = 1.;

    if (y < 0) {
        signy = -1;
        y = -y;
    }

    VecX<DIM> ret;
    double g, h;

    if (y >= my_gamma(DIM) * x.norm()) {
        //top
        g = gtop(x, y) / tau(DIM);
        h = htop(x, y);
    }
    else {
        //bot
        g = gbot(x, y);
        h = hbot(x, y) / rho(DIM);
    }
    for (int i = 0; i < DIM - 1; ++i) {
        ret[i] = x[i] * g;
    }
    ret[DIM - 1] = signy * h;

    return ret;

}

template<int DIM>
inline VecX<DIM> NCylinder2NCube(const VecX<DIM>& p) {
    if (p.norm() < 1E-16) {
        return p;
    }

    VecX<DIM> ret;
    VecX<DIM-1> x;
    for (int i = 0; i < DIM-1; ++i) {
        x[i] = p[i];
    }
    VecX<DIM - 1> newX = NBall2NCube(x);
    for (int i = 0; i < DIM-1; ++i) {
        ret[i] = newX[i];
    }
    ret[DIM - 1] = p[DIM - 1];
    return ret;
}

inline VecX<(int)1> NCylinder2NCube(const VecX<(int)1>& p){
    return p;
}

inline VecX<(int)2> NCylinder2NCube(const VecX<(int)2>& p){
    return p;
}

//Time taken by gsl: 110690 ms
//Time taken by nogsl: 17301 ms
template<int DIM>
inline VecX<DIM> NBall2NCube(const VecX<DIM>& p) {
    return NCylinder2NCube(NBall2NCylinder(p));
}
inline VecX<2> disk2square(const VecX<2>& p){
    return NBall2NCube(p);
}

inline VecX<1> NCylinder2NBall(const VecX<1>& p){
    return p;
}
template <int N>
inline VecX<N> NCylinder2NBall(const VecX<N>& p){
    int d = N;

    double y = p[d - 1];
    double signy = (y < 0. ? -1. : 1.);
    y = std::abs(y);
    double normx = std::sqrt(p.norm2() - y*y);

    VecX<N> res;
    if (y > normx){
        double xScale = y / normx;
        double invTau = inverseTau(tau(d) / xScale, d);
        double divider = std::sqrt(1. + invTau * invTau);

        for (int i = 0; i < d - 1; ++i){
            res[i] = p[i] * xScale / divider;
        }
        res[d-1] = signy * y * invTau / divider;
    } else {
        double invRho = inverseRho(rho(d) * y / normx, d);
        double divider = std::sqrt(1. + invRho * invRho);

        for (int i = 0; i < d - 1; ++i){
            res[i] = p[i] / divider;
        }
        res[d-1] = signy * normx * invRho / divider;
    }

    return res;
}

template <int N>
VecX<N> NCube2NBall(const VecX<N>& p);

inline VecX<1> NCube2NCylinder(const VecX<1>& p){
    return p;
}
inline VecX<2> NCube2NCylinder(const VecX<2>& p){
    return p;
}
template <int N>
inline VecX<N> NCube2NCylinder(const VecX<N>& p){
    VecX<N-1> x;
    for (int i = 0; i < N-1; ++i){
        x[i] = p[i];
    }
    x = NCube2NBall(x);
    VecX<N> ret;
    for (int i = 0; i < N-1; ++i){
        ret[i] = x[i];
    }
    ret[N - 1] = p[N - 1];
    return ret;
}


inline VecX<1> NCube2NBall(const VecX<1>& p){
    return p;
}
template <int N>
inline VecX<N> NCube2NBall(const VecX<N>& p){
    return NCylinder2NBall(NCube2NCylinder(p));
}

#endif //SLICEDOPTIM_MAPPING_H
