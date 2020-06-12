//
//

#ifndef SLICEDOPTIM_NBALLRADONMANAGER_H
#define SLICEDOPTIM_NBALLRADONMANAGER_H

#include <vector>
#include <string>
#include <iostream>


/**
 * Handles all information required for N dimensionnal ball Radon transform CDF inversion
 */
class NBallRadonManager {
private:
    int N;
    double normFactor;
    std::vector<double> coefs;
    std::vector<double> CDFcoefs;
    double CDFoffset;

    /**
     * Computes the factor by which the radonTransform must multiplied in order to integrate to 1
     *
     * @param d Dimension in which we are working
     * @return
     */
    double computeNormalisationFactor(int d);

public:
    NBallRadonManager(int dimension);

    /**
     * The Radon transform of the N Ball at position x is given by
     * \f$ R(x) = \frac{\pi^{d/2}}{\Gamma\bigl(\frac d 2 + 1\bigr)} \sqrt{1-x^2}^{d-1} \f$
     * @param x Position at which to compute the Radon Transform
     * @return The Radon transform of the N dimensional ball for position x along the direction
     */
    double radon(double x) const;

    /**
     * In odd dimensions the Radon Transform is a polynom which integrates to an other polynom.
     * In even dimensions it is more difficult to integrate.
     * That give us the following closed formula:
     * \f$C(x)=\begin{cases}
     * \sum\limits_{k=0}^{\frac{d-1}{2}}
     *    (-1)^k \left( \begin{smallmatrix} \tfrac{d-1}2\\ k \end{smallmatrix} \right)
     *		\frac{ x^{2k+1}}{2k+1} & \text{if $d$ is odd,}\\[4mm]
     * \frac{ \sqrt{\pi} \, \Gamma\bigl(\frac{1+d}{2} \bigr)}
     *   {2 \;\Gamma \bigl(1 + \frac{d}{2}\bigr)}
     *   +  {}_{\scriptscriptstyle 2\!}F_{\scriptscriptstyle1}\bigl(\frac{1}{2}, \frac{1-d}{2}, \frac{3}{2}, x^2\bigr) \, x & \text{if $d$ is even,}
     * \end{cases} \f$
     *
     * @param x Position at which to compute the Radon Transform's CDF
     * @return The CDF of the Radon transform of the N dimensional ball for position x along the direction
     */
    double radonCDF(double x) const;

    /**
     * Inverse the CDF of the radon transform of a N Ball using gradient descent.
     *
     * @param x Value for which to inverse the CDF
     * @return Value for which the Radon transform's CDF is x
     */
    double inverseCDF(double x) const;

};

#endif //SLICEDOPTIM_NBALLRADONMANAGER_H
