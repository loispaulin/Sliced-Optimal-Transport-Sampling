//
//

#ifndef SLICEDOPTIM_NBALLSLICEDOPTIMALTRANSPORT_H
#define SLICEDOPTIM_NBALLSLICEDOPTIMALTRANSPORT_H

#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <cstring>
#include <set>
#include <map>
#include "NBallRadonManager.h"
#include "../Math/VecX.h"
#include "../Math/myMath.h"
#include "../Tools/iopointset.h"
#include "../Tools/mapping.h"
#include "../Tools/my_utility.h"


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


/**
 * Choose \p m directions in N dimension N being defined by the dimention of the content of directions.
 *
 * @param directions Table of directions to output. Must be initialized with \p m VECTYPE of the disired dimension
 * @param seed Seed for the random generator. Only applied once
 * @param projectionsIndices Array containing the indices associated to the projections.
 * @param projectionsCumulatedSlicesNumber Array of the cumulative number of slices associated to the projections.
 */
template <class VECTYPE>
inline void newChooseDirectionsND(std::vector<VECTYPE>& directions, std::mt19937& generatorND, const std::vector<std::vector<uint>>& projectionsIndices, const std::vector<size_t>& projectionsCumulatedSlicesNumber){

    static std::normal_distribution<>normalND;

    const uint nbDims = directions.front().dim();
    const uint nbProjs = projectionsIndices.size();

    const uint nbDirs = directions.size();

    //Initializing directions to zero.
#pragma omp parallel for
    for(uint idDir = 0u; idDir < nbDirs; idDir++)
        std::memset(&directions[idDir][0], 0, nbDims * sizeof(directions[idDir][0]));

    //Building directions.
    for(uint idProj = 0u; idProj < nbProjs; ++idProj)
    {
//std::cerr << "proj " << idProj << std::endl;
        const uint firstInd = (idProj == 0u ? 0u : projectionsCumulatedSlicesNumber[idProj - 1u]);
        const uint lastInd = projectionsCumulatedSlicesNumber[idProj];
        const uint nbSlices = lastInd - firstInd;
        const uint effectiveNbDims = projectionsIndices[idProj].size();
//std::cerr << "firstInd " << firstInd << std::endl;
//std::cerr << "lastInd " << lastInd << std::endl;

        //uint i(0u);
        //for(auto indice : projectionsIndices[idProj])
        //    std::cerr << "indice " << i++ << ": " << indice << std::endl;

        const double pradius = 0.24 * std::pow(nsphereArea(1, effectiveNbDims) / (nbSlices * nballVolume(1, effectiveNbDims - 1)), (nbDims - 1) * 0.5);
        const double cosPradius = 1. / std::sqrt(1. + pradius * pradius);

//#pragma omp parallel for
        for(uint idDir = firstInd; idDir < lastInd; idDir++)
        {
            do
            {
                for(uint idProjDim(0u); idProjDim < effectiveNbDims; idProjDim++)
                    directions[idDir][projectionsIndices[idProj][idProjDim]] = normalND(generatorND);

                directions[idDir].normalize();
            }
            while(effectiveNbDims > 1 && !testPoissonND(directions[idDir], directions.begin() + firstInd, directions.begin() + idDir, cosPradius));
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
 * To be correct Projection wise, points are given in a cube domain and transport will be computed on ball mapping
 * of targeted projection
 *
 * @param dir Direction \f$ \theta \f$
 * @param projectionIndices Table of the dimensions of the targeted projection
 * @param points Table containing the points \f$ x_j \f$ in cube domain
 * @param pos Table containing the optimal solution in 1D
 * @param shift Output the 1D shift to apply to minimize transport cost
 * @param pointsProject Memory buffer used to store the 1D projection of points. Must be the same size as \p pointsi
 * @param cylinderPoints
 */
template<class VECTYPE>
inline void newSlicedStepNBall(
        const VECTYPE& dir,
        const std::vector<uint>& projectionIndices,
        const std::vector<VECTYPE>& points,
        const std::vector<std::vector<double>>& pos,
        std::vector<VECTYPE>& shift,
        std::vector<std::pair<double, int>>& pointsProject,
        std::vector<VecXDynamic> cylinderPoints
) {
    //Number of dimensions in total
    const uint nbDims = dir.dim();
    //Number of dimensions used in current projection
    const uint effectiveNbDims = projectionIndices.size();
    //Number of points.
    const size_t nbPoints = points.size();

    //Compute direction vector restricted to these dimensions
    VecXDynamic effectiveDir(effectiveNbDims);
    for(uint i(0u); i < effectiveNbDims; i++)
        effectiveDir[i] = dir[projectionIndices[i]];

   //Computes points projection and sort them along given direction
    for(size_t idPoint = 0u; idPoint < nbPoints; ++idPoint)
    {

        //Restrict cube to used dimensions
        VecXDynamic effectiveP(effectiveNbDims);
        for(uint idProjDim(0u); idProjDim < effectiveNbDims; idProjDim++)
            effectiveP[idProjDim] = points[idPoint][projectionIndices[idProjDim]];

        //Cube to Cylinder == (projection of cube) to ball
        cylinderPoints[idPoint] = NCube2NBall(effectiveP);

        //Compute projection in cylindrical domain
        pointsProject[idPoint].first = effectiveDir * cylinderPoints[idPoint];
        pointsProject[idPoint].second = int(idPoint);
    }

    //Do 1D transport in cylindrical domain
    std::sort(pointsProject.begin(), pointsProject.end());

    //Computes required shift to optimize 1D optimal transport
    for(size_t idPoint = 0u; idPoint < nbPoints; ++idPoint)
    {
        //Compute shifting
        const double s = pos[effectiveNbDims - 1][idPoint] - pointsProject[idPoint].first;

        const int indice = pointsProject[idPoint].second;

        //Apply shift in cylindrical domain
        cylinderPoints[indice] += s * effectiveDir;
        //Compute projected cube coordinates
        VecXDynamic effectiveP = NBall2NCube(cylinderPoints[indice]);

        //Get shift in cube coordinates
        VECTYPE& pointShift = shift[indice];
        //Initializing whole cube vector to zero before setting projection shift components
        std::memset(&pointShift[0], 0, nbDims * sizeof(pointShift[0]));
        for(uint idProjDim(0u); idProjDim < effectiveNbDims; idProjDim++)
        {
            const uint idDim = projectionIndices[idProjDim];
            pointShift[idDim] = effectiveP[idProjDim] - points[indice][idDim];
        }
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
 * @param projectionsIndices
 * @param projectionsSizesCylinders
 * @param projectionsCumulatedSlicesNumber
 * @param pos Table containing the 1D positions optimizing transport in 1D
 * @param shifts Used to avoid having to allocate uge chunks of memory. Shifts for each direction and each point.
 * @param finalShift Used to avoid having to allocate huge chunks of memory Must be a vector of same size as \p pointsOut.
 * @param pointsProject Used to avoid having to allocate uge chunks of memory. Must be a vector of size m containing vectors of same size as \p pointsOut.
 * @param cylindersPoints Used to avoid having to allocate huge chunks of memory.
 * @param leastSquareMethod Enables Least square advection method
 * @param nesterov Enables Nesterov gradient descent method
 * @param lambdas
 * @param yS
 * @param displacement
 */
template <class VECTYPE>
inline void     newSlicedOptimalTransportBatch(std::vector<VECTYPE>& pointsOut,
                                               const std::vector<VECTYPE>& directions,
                                               const std::vector<std::vector<uint>>& projectionsIndices,
                                               const std::map<uint, uint>& projectionsCylinders,
                                               const std::vector<size_t>& projectionsCumulatedSlicesNumber,
                                               const std::vector<double>& projectionsWeights,
                                               const std::vector<std::vector<double>>& pos,
                                               std::vector<std::vector<VECTYPE>>& shifts,
                                               std::vector<VECTYPE>& finalShift,
                                               std::vector<std::vector<std::pair<double, int>>>& pointsProject,
                                               std::vector<std::vector<VecXDynamic>>& cylindersPoints,
                                               std::vector<double>& displacement){
    const size_t nbPoints = pointsOut.size();
    const uint nbDims = directions.front().dim();
    const uint nbProjs = projectionsIndices.size();

    //Compute the shift along each direction
    for(uint idProj(0u); idProj < nbProjs; idProj++)
    {
        const uint firstInd = (idProj == 0u ? 0u : projectionsCumulatedSlicesNumber[idProj - 1u]);
        const uint lastInd = projectionsCumulatedSlicesNumber[idProj];

        const uint cylinderId = projectionsCylinders.at(idProj);

#pragma omp parallel for
        for(uint idDir = firstInd; idDir < lastInd; ++idDir)
            newSlicedStepNBall(directions[idDir], projectionsIndices[idProj], pointsOut, pos, shifts[idDir], pointsProject[idDir], cylindersPoints[cylinderId]);
    }

    //Accumulate shift from all directions
#pragma omp parallel for
    for(size_t idPoint = 0u; idPoint < nbPoints; ++idPoint)
    {
        VECTYPE& shift(finalShift[idPoint]);
        std::memset(&shift[0], 0, nbDims * sizeof(shift[0]));

        for(uint idProj(0u); idProj < nbProjs; idProj++)
        {
            const uint nbProjDims = projectionsIndices[idProj].size();
            const uint firstInd = (idProj == 0u ? 0u : projectionsCumulatedSlicesNumber[idProj - 1u]);
            const uint lastInd = projectionsCumulatedSlicesNumber[idProj];
            const uint nbDirProj = lastInd - firstInd;

            for(uint idDir(firstInd); idDir < lastInd; idDir++)
                for(uint idProjDim(0u); idProjDim < nbProjDims; idProjDim++)
                {
                    const uint idDim = projectionsIndices[idProj][idProjDim];
                    shift[idDim] += shifts[idDir][idPoint][idDim] * projectionsWeights[idProj] / nbDirProj;
                }
        }

    }

    //Compute next step
    //Displace points according to accumulated shift
#pragma omp parallel for
    for(size_t idPoint = 0u; idPoint < nbPoints; ++idPoint) {
        pointsOut[idPoint] += finalShift[idPoint];
        displacement[idPoint] = finalShift[idPoint].norm();
    }
}


/**
 *  \brief Computes an optimized point set to uniformly sample the [-1, 1]^d Cube using sliced optimal transport
 *
 * @param pointsIn Contains input ND points
 * @param pointsOut Contains optimized ND points
 * @param nbIter Number of iterations
 * @param nbSlices Total number of slices per iteration
 * @param seed random seed
 * @param projectionsMasks
 * @param projectionsSlicesNumber
 * @param leastSquareMethod
 */
template <class VECTYPE>
inline void newSlicedOptimalTransport(const std::vector<VECTYPE>& pointsIn,
                                      std::vector<VECTYPE>& pointsOut,
                                      int nbIter,
                                      int nbSlices,
                                      std::mt19937& generatorND,
                                      const std::vector<std::vector<bool> >& projectionsMasks,
                                      const std::vector<size_t>& projectionsCumulatedSlicesNumber,
                                      const std::vector<double>& projectionsWeights){

    const uint nbDims = pointsIn.front().dim();
    const size_t nbSamples = pointsIn.size();
    const uint nbProjs = projectionsMasks.size();

    pointsOut = pointsIn;

    //Compute optimal 1D positions for given number of points and each possible size of projection
    std::vector<std::vector<double>> pos(nbDims, std::vector<double>(nbSamples));
    for(uint indDim = 0u; indDim < nbDims; ++indDim)
        getInverseRadonNBall(int(indDim) + 1, nbSamples, pos[indDim]);

//std::cerr << "inverse done" << std::endl;
//std::cerr << "projections number " << projectionsMasks.size() << std::endl;

    //Build indices associated to projections.
    std::vector<std::vector<uint>> projectionsIndices(projectionsMasks.size(), std::vector<uint>());
    for(uint indProj = 0u; indProj < nbProjs; indProj++)
        for(uint indDim = 0u; indDim < nbDims; indDim++)
            if(projectionsMasks[indProj][indDim])
                projectionsIndices[indProj].emplace_back(indDim);

    //Storage for sorting information of the slice 1D projected points.
    std::vector<std::vector<std::pair<double, int>>> pointsProject(nbSlices, std::vector<std::pair<double, int>>(nbSamples));
    //Accumulation shift to be applied later on
    std::vector<std::vector<VECTYPE>> shifts(nbSlices, std::vector<VECTYPE>(nbSamples, VECTYPE(nbDims)));
    std::vector<VECTYPE> finalShift(nbSamples, VECTYPE(nbDims));

    //Set of projection sizes.
    std::set<uint> projectionsSizesSet;
    for(uint idProj = 0u; idProj < nbProjs; idProj++)
        projectionsSizesSet.insert(projectionsIndices[idProj].size());
    const uint nbProjsSizes = projectionsSizesSet.size();

    //std::cerr << "Projections sizes set:" << std::endl;
    //for(const uint size : projectionsSizesSet)
    //    std::cerr << size << std::endl;

    //Storage for cylinders projected points.
    std::vector<std::vector<VecXDynamic>> cylindersPoints(nbProjsSizes);
    uint idProjSize(0u);
    for(const uint size : projectionsSizesSet)
    {
        cylindersPoints[idProjSize++] = std::vector<VecXDynamic>(nbSamples, VecXDynamic(size));
    }

    //Map of projection--cylinder associations.
    std::map<uint, uint> projectionsCylinders;
    for(uint idProj = 0u; idProj < nbProjs; idProj++)
    {
        const uint projSize = projectionsIndices[idProj].size();
        auto it = std::find(projectionsSizesSet.begin(), projectionsSizesSet.end(), projSize);
        const uint cylinderId = std::distance(projectionsSizesSet.begin(), it);
        //std::cerr << "Projection " << idProj << " uses " << projSize << " dimensions and will use the cylinder " << cylinderId << std::endl;
        projectionsCylinders.emplace(std::make_pair(idProj, cylinderId));
    }

    std::vector<double> displacement(nbSamples);
    //Iterate 1D Optimal Transport
    std::vector<VECTYPE> directions(nbSlices, VECTYPE(nbDims));
    for(uint i(0u); i < static_cast<uint>(nbIter); i++){
        newChooseDirectionsND(directions, generatorND, projectionsIndices, projectionsCumulatedSlicesNumber);
        newSlicedOptimalTransportBatch(pointsOut, directions, projectionsIndices, projectionsCylinders,
                                       projectionsCumulatedSlicesNumber, projectionsWeights, pos, shifts, finalShift,
                                       pointsProject, cylindersPoints, displacement);
    }

}


template <class VECTYPE>
inline bool testPoissonND(const VECTYPE& v, typename std::vector<VECTYPE>::const_iterator begin, typename std::vector<VECTYPE>::const_iterator end, double radius){
    bool test = true;
    for (typename std::vector<VECTYPE>::const_iterator it = begin; it != end && test; ++it){
        test = (v * (*it) < radius);
    }
    return test;
}


#endif //SLICEDOPTIM_NBALLSLICEDOPTIMALTRANSPORT_H
