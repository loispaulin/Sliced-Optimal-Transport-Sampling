//
//

#ifndef SLICEDOPTIM_MY_UTILITY_H
#define SLICEDOPTIM_MY_UTILITY_H

#include <random>
#include <functional>
#include <string>
#include "../Math/VecX.h"

template<class VECTYPE, class RandomGenerator>
inline VECTYPE randomVectorInBall(int dim, RandomGenerator &engine) {
    std::normal_distribution<double> normalDistribution(0, 1);
    std::uniform_real_distribution<double> unif(0, 1);
    VECTYPE v(dim);

    for (int j = 0; j < dim; ++j) {
        v[j] = normalDistribution(engine);
    }
    v.normalize();
    v *= std::pow(unif(engine), 1. / dim);

    return v;
}

double clamp(double v, double min, double max);

double inverseFunction(std::function<double(double)> &f, std::function<double(double)> &df, double v);

#endif //SLICEDOPTIM_MY_UTILITY_H
