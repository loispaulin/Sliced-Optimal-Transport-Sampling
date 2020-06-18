//
//

#include <cmath>
#include "../Math/myMath.h"
#include "../Tools/my_utility.h"
#include "mapping.h"

/**
 * Returns the gamma value used in the mapping for given dimention \p dim
 * @param dim
 * @return Mapping gamma value
 */
double my_gamma(int dim){

    static double gammas[65] = {0, 0, 1, 0.8944271909999136, 0.8213530892079441, 0.7666031370294691, 0.7234249021341935, 0.6881297272460598, 0.6585046305043634, 0.6331279880529028, 0.6110376442187433, 0.5915550708605355, 0.5741855699257183, 0.5585584223422113, 0.5443893320612817, 0.5314559515561526, 0.5195814118815603, 0.5086229291980788, 0.4984637317388613, 0.4890072179012491, 0.4801726494508217, 0.4718919233755111, 0.4641071160092422, 0.4567685894686779, 0.4498335138119681, 0.4432647008230895, 0.4370296743575151, 0.4310999223518928, 0.4254502898314481, 0.4200584824268004, 0.4149046572977357, 0.409971083769389, 0.4052418600097225, 0.4007026750911621, 0.3963406080484276, 0.3921439572921156, 0.3881020950921822, 0.3842053428689888, 0.380444863828941, 0.3768125702132689, 0.3733010427315929, 0.3699034605, 0.3666135399, 0.3634254805, 0.360333918, 0.3573338828, 0.3544207626, 0.35159027, 0.3488384132, 0.3461614704, 0.3435559656, 0.3410186493, 0.3385464788, 0.3361366006, 0.3337863383, 0.3314931739, 0.3292547418, 0.3270688145, 0.3249332784, 0.32284617, 0.3208055833, 0.3188097766, 0.3168570517, 0.3149458485, 0.3130746005};

    if (dim < 2 || dim > 64){
        std::cerr << "Ball <-> Cylinder mapping not yet implemented in dimension " << dim << std::endl;
        exit(1);
    } else {
        return gammas[dim];
    }
}

/**
 * Returns tau_dim(gamma) for mapping usage
 * @param dim
 * @return
 */
double tau(int dim){

    static double taus[65] = {0, 0, 0.7853981633974483, 0.8164965809277261, 0.8382695966098712, 0.8545740127924695, 0.8673491949880974, 0.8776916965664365, 0.8862745508336506, 0.8935367660649958, 0.899778490075909, 0.9052127722373648, 0.9099955397200091, 0.9142438593879563, 0.9180475385727878, 0.9214767582425344, 0.9245872485188505, 0.9274238899409514, 0.9300232766197379, 0.9324155772091874, 0.9346259101418143, 0.9366753760830803, 0.9385818441258782, 0.9403605581967132, 0.9420246102751126, 0.9435853136309325, 0.9450525000901884, 0.9464347589278826, 0.9477396304476875, 0.9489737640514276, 0.9501430482355193, 0.9512527182118499, 0.9523074455581219, 0.9533114133316983, 0.9542683793486204, 0.9551817297641813, 0.9560545246599033, 0.9568895370081494, 0.9576892861133089, 0.9584560664420457, 0.9591919725533894, 0.9598989208000001, 0.960578668, 0.9612328283, 0.9618628873, 0.9624702149, 0.9630560765, 0.9636216432, 0.964168, 0.9646961541, 0.965207042, 0.9657015353, 0.9661804463, 0.9666445337, 0.9670945059, 0.9675310264, 0.9679547163, 0.9683661587, 0.9687658998999999, 0.9691544549, 0.9695323071, 0.9698999117, 0.9702577002, 0.9706060771, 0.970945422};

    if (dim < 2 || dim > 64){
        std::cerr << "Ball <-> Cylinder mapping not yet implemented in dimension " << dim << std::endl;
        exit(1);
    } else {
        return taus[dim];
    }

}

/**
 * Returns rho_dim(gamma) for mapping usage
 * @param dim
 * @return
 */
double rho(int dim){

    static double rhos[65] = {0, 0, 0.7853981633974483, 0.6666666666666666, 0.5890486225480868, 0.5333333333333324, 0.4908738521234047, 0.4571428571428578, 0.4295146206079796, 0.4063492063492068, 0.3865631585471811, 0.369408369408369, 0.3543495620015831, 0.3409923409923414, 0.3290388790014698, 0.3182595182595185, 0.3084739490638782, 0.2995383701266054, 0.2913365074492179, 0.2837731927515205, 0.2767696820767574, 0.270260183572876, 0.264189241982359, 0.2585097408088374, 0.2531813568997612, 0.2481693511764841, 0.2434436124036166, 0.2389778937254978, 0.2347491976749156, 0.2307372767004782, 0.2269242244190854, 0.2232941387424428, 0.2198328424059889, 0.216527649689546, 0.2133671705705185, 0.210341145412814, 0.2074403047213377, 0.204656249590948, 0.2019813493339337, 0.1994086534480554, 0.1969318156005856, 0.1945450278, 0.1922429628, 0.1900207248, 0.1878738046, 0.185798042, 0.1837895915, 0.1818448922, 0.1799606416, 0.1781337719, 0.1763614288, 0.1746409529, 0.1729698629, 0.1713458405, 0.1697667173, 0.1682304616, 0.1667351687, 0.1652790502, 0.1638604244, 0.1624777107, 0.1611294174, 0.1598141414, 0.1585305558, 0.1572774097, 0.1560535159};

    if (dim < 2 || dim > 64){
        std::cerr << "Ball <-> Cylinder mapping not yet implemented in dimension " << dim << std::endl;
        exit(1);
    } else {
        return rhos[dim];
    }

}

/**
 * Maps the N-Ball into a N-Cylinder
 * @param p
 * @return
 */


double tau(double lambda, int dim){
    //std::cout << integrateSinPower(0, std::atan(1. / lambda), dim - 2) << std::endl;
    //std::cout << "{tau[" << lambda << ", " << dim << "], " << std::pow( (dim-1) * integrateSinPower(0, std::atan(1. / lambda), dim - 2), 1. / (dim-1)) << "}," << std::endl;
    return std::pow( (dim-1) * integrateSinPower(0, std::atan(1. / lambda), dim - 2), 1. / (dim-1));
}
double dtau(double lambda, int dim){
    return -(1. / (std::pow(tau(lambda, dim), dim - 2) * std::pow(1. + lambda * lambda, dim * 0.5)));
}
double inverseTau(double v, int dim){

    std::function<double(double)> t = [&dim](double l) { return tau(l, dim); };
    std::function<double(double)> dt = [&dim](double l) { return dtau(l, dim); };

    return inverseFunction(t, dt, v);

}

double rho(double lambda, int dim){
    //std::cout << "{rho[" << lambda << ", " << dim << "], " << integrateCosPower(0, std::atan(lambda), dim - 2) << "}," << std::endl;
    return integrateCosPower(0, std::atan(lambda), dim - 2);
}
double drho(double lambda, int dim){
    return 1. / std::pow(1 + lambda * lambda, dim * 0.5);
}
double inverseRho(double v, int dim){

    std::function<double(double)> r = [&dim](double l) { return rho(l, dim); };
    std::function<double(double)> dr = [&dim](double l) { return drho(l, dim); };

    return inverseFunction(r, dr, v);

}


double leftgs(double gamma, int dim){
    return (dim - 1) * integrateSinPower(0, std::atan(1. / gamma), dim - 2);
}
double rightgs(double gamma, int dim){
    return integrateCosPower(0, std::atan(gamma), dim - 2);
}
double evalgs(double gamma, int dim){
    return leftgs(gamma, dim) - rightgs(gamma, dim);
}

double solveForGamma(int dim){

    double eps = std::pow(10., -14);

    double mini = 0.01;
    double maxi = 0.99;

    while (maxi - mini > eps){
        double mid = 0.5 * (maxi + mini);
        double v = evalgs(mid, dim);
        if (v < 0){
            maxi = mid;
        } else {
            mini = mid;
        }
    }

    return 0.5 * (maxi + mini);
}
