//
//

#ifndef SLICEDOPTIM_NBALLSLICEDOPTIMALTRANSPORT_H
#define SLICEDOPTIM_NBALLSLICEDOPTIMALTRANSPORT_H

#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <cstring>
#include "NBallRadonManager.h"
#include "../Math/VecX.h"
#include "../Math/myMath.h"
#include "../Tools/iopointset.h"
#include "../Tools/mapping.h"
#include "../Tools/my_utility.h"

template <class VECTYPE>
inline bool testPoissonND(const VECTYPE& v, typename std::vector<VECTYPE>::const_iterator begin, typename std::vector<VECTYPE>::const_iterator end, double radius){
    bool test = true;
    for (typename std::vector<VECTYPE>::const_iterator it = begin; it != end && test; ++it){
        test = (v * (*it) < radius);
    }
    return test;
}

/**
 * Choose \p m directions in N dimension N being defined by the dimention of the content of directions.
 * Two selection methods are available. Either the direction are uniformly selected in ND
 * or one can force half of the them to lie in 2D planes to optimize the projection repartition as well.
 *
 * @param directions Table of directions to output. Must be initialized with \p m VECTYPE of the disired dimension
 * @param m Number of directions to pick
 * @param seed Seed for the random generator. Only applied once
 * @param projective If true half the directions will lie in 2D planes.
 */
template <class VECTYPE>
inline void chooseDirectionsND(std::vector<VECTYPE>& directions, int m, int seed){

    static std::mt19937 generatorND(seed);
    static std::normal_distribution<>normalND;

    int dim = directions.front().dim();

    double pradius = std::cos(0.99 * std::pow(nsphereArea(1, dim) / (m * nballVolume(1, dim - 1)), 1. / (dim - 1)));
    for (int k = 0; k < m; ++k){
        for (int j = 0; j < dim; ++j){
            directions[k][j] = normalND(generatorND);
        }
        directions[k].normalize();
        if (!testPoissonND(directions[k], directions.begin(), directions.begin() + k, pradius)){
            k -= 1;
        }
    }
}

/**
 * Compute optimal transport in 1D for direction \f$ \theta \f$ and \f$d_{j}\f$ being the 1D displacement of \f$\x^j\f$
 * that minimize the 1D sliced optimal transport along \f$ \theta \f$.
 *
 * Denoting $\sigma$ the permutations of the indices \f$\{j\}_{j=1..N}\f$ such that
 * \f$\bigl(\x^{\sigma(j)} \cdot \theta \bigr)_j\f$, is a sorted sequence of increasing values,
 * one can compute \f$d_{j}\f$ via
 * \f$ d_{j} = C_{\theta}^{-1}\left(\frac{\sigma(j)-\frac12}{N}\right)\,. \vspace*{-1mm}\f$
 *
 * @param dir Direction \f$ \theta \f$
 * @param points Table containing the points \f$ x_j \f$
 * @param pos Table containing the optimal solution in 1D
 * @param shift Output the 1D shift to apply to minimize transport cost
 * @param pointsProject Memory buffer used to store the 1D projection of points. Must be the same size as \p points
 */
template<class VECTYPE>
inline void slicedStepNBall(const VECTYPE& dir,
                            const std::vector<VECTYPE>& points,
                            const std::vector<double>& pos,
                            std::vector<double>& shift,
                            std::vector<std::pair<double, int>>& pointsProject)
{

    //Computes points projection and sort them along given direction
    for (size_t i = 0; i < pointsProject.size(); ++i){
        pointsProject[i].first = dir * points[i];
        pointsProject[i].second = int(i);
    }

    std::sort(pointsProject.begin(), pointsProject.end());

    //Computes required shift to optimize 1D optimal transport
    for (size_t i = 0; i < pointsProject.size(); ++i) {
        //Compute shifting
        double s = pos[i] - pointsProject[i].first;

        shift[pointsProject[i].second] = s;
    }

}

/**
 * Compute optimal transport in 1D for the \p directions and displace \f$ x_j \f$ by
 * \f$\pmb{\delta}^j \EqDef \frac{1}{K} \sum_{i=1}^K d_{i,j}\, \theta_i \vspace*{-1mm}\f$ with
 * with \f$ d_{i,j} \f$ being the displacement of \f$\x^j\f$ that minimize the 1D sliced optimal
 * transport along direction \f$ \theta_i \f$
 *
 * @param pointsOut Table containing the points \f$ x_j \f$
 * @param directions Table containing the \f$\theta_i\f$
 * @param pos Table containing the 1D positions optimizing transport in 1D
 * @param shift Used to avoid having to allocate uge chunks of memory. Must be a vector of size m containing vectors of same size as \p pointsOut.
 * @param finalShift Used to avoid having to allocate huge chunks of memory Must be a vector of same size as \p pointsOut.
 * @param pointsProject Used to avoid having to allocate uge chunks of memory. Must be a vector of size m containing vectors of same size as \p pointsOut.
 * @return the Wasserstein cost of the current iteration.
 */
template <class VECTYPE>
inline void slicedOptimalTransportBatch(std::vector<VECTYPE>& pointsOut,
                                 const std::vector<VECTYPE>& directions,
                                 const std::vector<double>& pos,
                                 std::vector<std::vector<double>>& shift,
                                 std::vector<VECTYPE>& finalShift,
                                 std::vector<std::vector<std::pair<double, int>>>& pointsProject)
 {

    int m = directions.size();
    int nbPoints = pointsOut.size();

    //Compute the shift along each direction
#pragma omp parallel for shared(directions, shift)
    for (int k = 0; k < m; ++k){
        for(double& v : shift[k]){
			v = 0.;
        }
        const VECTYPE& dir = directions[k];

        slicedStepNBall(dir, pointsOut, pos, shift[k], pointsProject[k]);
    }

        //Accumulate shift from all directions
#pragma omp parallel for
    for (int i = 0; i < nbPoints; ++i) {
        VECTYPE sh(finalShift[i].dim());
        memset(&sh[0], 0, finalShift[i].dim() * sizeof(sh[0]));
        for (int k = 0; k < m; ++k) {
            sh += shift[k][i] * directions[k];
        }
        finalShift[i] = sh;
        finalShift[i] /= m;
    }

    //Displace points according to accumulated shift
#pragma omp parallel for
    for (int i = 0; i < nbPoints; ++i) {
        pointsOut[i] += finalShift[i];
    }

}

/**
 *  \brief Get Optimal position at which to place \p nbSamples 1D sample to minimize OT cost to the Radon transform of a \p D dimensional ball
 *
 *  A common result in optimal transport is that 1D placement comes from the inverse of the CDF of target distribution
 * @param D Dimension of the Ball
 * @param nbSamples Number of samples to compute position for
 * @param pos Output buffer in which to put the positions
 */
inline void getInverseRadonNBall(int D, int nbSamples, std::vector<double>& pos){

    pos.resize(nbSamples);
    NBallRadonManager nbrm(D);

#pragma omp parallel for shared(pos)
    for (int i = 0; i < nbSamples; ++i){
        double p = (2. * i + 1.) / (2. * nbSamples);
        pos[i] = nbrm.inverseCDF(p);
    }
}

/**
 *  \brief Computes an optimized point set to uniformly sample the unit N-Ball using sliced optimal transport
 *
 * @param pointsIn Contains input ND points
 * @param pointsOut Contains optimized ND points
 * @param nbIter Number of iterations
 * @param m Number of slice per iteration
 * @param seed random seed
 * @return the Sliced Wasserstein distance between the samples and the uniform distribution.
 */

template <class VECTYPE>
inline void slicedOptimalTransportNBall(const std::vector<VECTYPE>& pointsIn,
                                        std::vector<VECTYPE>& pointsOut,
                                        int nbIter,
                                        int m,
                                        int seed)
{

    int N = pointsIn.front().dim();
    pointsOut = pointsIn;

    //Compute optimal 1D position for given number of points;
    int nbSamples = int(pointsOut.size());
    std::vector<double> pos(pointsOut.size());
    getInverseRadonNBall(N, nbSamples, pos);

    //Allocate Memory to compute projections in
    std::vector<std::vector<std::pair<double, int>>> pointsProject(m, std::vector<std::pair<double, int>>(pointsOut.size()));

    //Accumulation shift to be applied later on
    std::vector<std::vector<double>> shift(m, std::vector<double>(pointsOut.size()));
    std::vector<VECTYPE> finalShift(pointsOut.size(), VECTYPE(N));

    std::vector<VECTYPE> directions(m, VECTYPE(N));

    //Iterate 1D Optimal Transport
    for (int i = 0; i < nbIter; i += 1){
        chooseDirectionsND(directions, m, seed);

        slicedOptimalTransportBatch(pointsOut, directions, pos, shift, finalShift, pointsProject);
    }

}

/**
 * Compute the 1D transport cost of the from given ND \param points
 * to the distribution whose optimal 1D positions along dir are given in \param pos
 *
 * @tparam VECTYPE any vector type with * being the dot product and [] granting access to the coordinates value
 * @param points ND points
 * @param dir Projection direction
 * @param pos Sorted optimal 1D positions for target distribution along dir
 * @param pointsProject Pre allocated buffer used to sort the points's projection
 */
template <class VECTYPE>
inline double computeSlicedOTCost(const std::vector<VECTYPE>& points,
                                  const VECTYPE& dir,
                                  const std::vector<double>& pos,
                                  std::vector<double>& pointsProject
                                  ){
    double cost = 0.;

    for (size_t i = 0; i < pointsProject.size(); ++i){
        pointsProject[i] = points[i] * dir;
    }
    sort(pointsProject.begin(), pointsProject.end());

    //Place them at optimal places
    for (size_t i = 0; i < pointsProject.size(); ++i) {
        cost += (pointsProject[i] - pos[i]) * (pointsProject[i] - pos[i]);
    }
    return cost / points.size();
}

#endif //SLICEDOPTIM_NBALLSLICEDOPTIMALTRANSPORT_H
