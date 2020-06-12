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
    switch (dim){case 2:
            return (1.);
        case 3:
            return 0.8944271909999136;
        case 4:
            return 0.8213530892079441;
        case 5:
            return 0.7666031370294691;
        case 6:
            return 0.7234249021341935;
        case 7:
            return 0.6881297272460598;
        case 8:
            return 0.6585046305043634;
        case 9:
            return 0.6331279880529028;
        case 10:
            return 0.6110376442187433;
        case 11:
            return 0.5915550708605355;
        case 12:
            return 0.5741855699257183;
        case 13:
            return 0.5585584223422113;
        case 14:
            return 0.5443893320612817;
        case 15:
            return 0.5314559515561526;
        case 16:
            return 0.5195814118815603;
        case 17:
            return 0.5086229291980788;
        case 18:
            return 0.4984637317388613;
        case 19:
            return 0.4890072179012491;
        case 20:
            return 0.4801726494508217;
        case 21:
            return 0.4718919233755111;
        case 22:
            return 0.4641071160092422;
        case 23:
            return 0.4567685894686779;
        case 24:
            return 0.4498335138119681;
        case 25:
            return 0.4432647008230895;
        case 26:
            return 0.4370296743575151;
        case 27:
            return 0.4310999223518928;
        case 28:
            return 0.4254502898314481;
        case 29:
            return 0.4200584824268004;
        case 30:
            return 0.4149046572977357;
        case 31:
            return 0.409971083769389;
        case 32:
            return 0.4052418600097225;
        case 33:
            return 0.4007026750911621;
        case 34:
            return 0.3963406080484276;
        case 35:
            return 0.3921439572921156;
        case 36:
            return 0.3881020950921822;
        case 37:
            return 0.3842053428689888;
        case 38:
            return 0.380444863828941;
        case 39:
            return 0.3768125702132689;
        case 40:
            return 0.3733010427315929;
        case 41:
            return 0.3699034605;
        case 42:
            return 0.3666135399;
        case 43:
            return 0.3634254805;
        case 44:
            return 0.360333918;
        case 45:
            return 0.3573338828;
        case 46:
            return 0.3544207626;
        case 47:
            return 0.35159027;
        case 48:
            return 0.3488384132;
        case 49:
            return 0.3461614704;
        case 50:
            return 0.3435559656;
        case 51:
            return 0.3410186493;
        case 52:
            return 0.3385464788;
        case 53:
            return 0.3361366006;
        case 54:
            return 0.3337863383;
        case 55:
            return 0.3314931739;
        case 56:
            return 0.3292547418;
        case 57:
            return 0.3270688145;
        case 58:
            return 0.3249332784;
        case 59:
            return 0.32284617;
        case 60:
            return 0.3208055833;
        case 61:
            return 0.3188097766;
        case 62:
            return 0.3168570517;
        case 63:
            return 0.3149458485;
        case 64:
            return 0.3130746005;
        default:
            std::cerr << "Ball <-> Cylinder mapping not yet implemented in dimension " << dim << std::endl;
            exit(1);
    }
}

/**
 * Returns tau_dim(gamma) for mapping usage
 * @param dim
 * @return
 */
double tau(int dim){
    switch (dim){case 2:
            return (M_PI_4);
        case 3:
            return (M_SQRT2 / sqrt(3));
        case 4:
            return 0.8382695966098712;
        case 5:
            return 0.8545740127924695;
        case 6:
            return 0.8673491949880974;
        case 7:
            return 0.8776916965664365;
        case 8:
            return 0.8862745508336506;
        case 9:
            return 0.8935367660649958;
        case 10:
            return 0.899778490075909;
        case 11:
            return 0.9052127722373648;
        case 12:
            return 0.9099955397200091;
        case 13:
            return 0.9142438593879563;
        case 14:
            return 0.9180475385727878;
        case 15:
            return 0.9214767582425344;
        case 16:
            return 0.9245872485188505;
        case 17:
            return 0.9274238899409514;
        case 18:
            return 0.9300232766197379;
        case 19:
            return 0.9324155772091874;
        case 20:
            return 0.9346259101418143;
        case 21:
            return 0.9366753760830803;
        case 22:
            return 0.9385818441258782;
        case 23:
            return 0.9403605581967132;
        case 24:
            return 0.9420246102751126;
        case 25:
            return 0.9435853136309325;
        case 26:
            return 0.9450525000901884;
        case 27:
            return 0.9464347589278826;
        case 28:
            return 0.9477396304476875;
        case 29:
            return 0.9489737640514276;
        case 30:
            return 0.9501430482355193;
        case 31:
            return 0.9512527182118499;
        case 32:
            return 0.9523074455581219;
        case 33:
            return 0.9533114133316983;
        case 34:
            return 0.9542683793486204;
        case 35:
            return 0.9551817297641813;
        case 36:
            return 0.9560545246599033;
        case 37:
            return 0.9568895370081494;
        case 38:
            return 0.9576892861133089;
        case 39:
            return 0.9584560664420457;
        case 40:
            return 0.9591919725533894;
        case 41:
            return 0.9598989208;
        case 42:
            return 0.960578668;
        case 43:
            return 0.9612328283;
        case 44:
            return 0.9618628873;
        case 45:
            return 0.9624702149;
        case 46:
            return 0.9630560765;
        case 47:
            return 0.9636216432;
        case 48:
            return 0.964168;
        case 49:
            return 0.9646961541;
        case 50:
            return 0.965207042;
        case 51:
            return 0.9657015353;
        case 52:
            return 0.9661804463;
        case 53:
            return 0.9666445337;
        case 54:
            return 0.9670945059;
        case 55:
            return 0.9675310264;
        case 56:
            return 0.9679547163;
        case 57:
            return 0.9683661587;
        case 58:
            return 0.9687658999;
        case 59:
            return 0.9691544549;
        case 60:
            return 0.9695323071;
        case 61:
            return 0.9698999117;
        case 62:
            return 0.9702577002;
        case 63:
            return 0.9706060771;
        case 64:
            return 0.970945422;
        default:
            std::cerr << "Ball <-> Cylinder mapping not yet implemented in dimension " << dim << std::endl;
            exit(1);
    }
}

/**
 * Returns rho_dim(gamma) for mapping usage
 * @param dim
 * @return
 */
double rho(int dim){
    switch (dim){
        case 2:
            return (M_PI_4);
        case 3:
            return (2./3.);
        case 4:
            return 0.5890486225480868;
        case 5:
            return 0.5333333333333324;
        case 6:
            return 0.4908738521234047;
        case 7:
            return 0.4571428571428578;
        case 8:
            return 0.4295146206079796;
        case 9:
            return 0.4063492063492068;
        case 10:
            return 0.3865631585471811;
        case 11:
            return 0.369408369408369;
        case 12:
            return 0.3543495620015831;
        case 13:
            return 0.3409923409923414;
        case 14:
            return 0.3290388790014698;
        case 15:
            return 0.3182595182595185;
        case 16:
            return 0.3084739490638782;
        case 17:
            return 0.2995383701266054;
        case 18:
            return 0.2913365074492179;
        case 19:
            return 0.2837731927515205;
        case 20:
            return 0.2767696820767574;
        case 21:
            return 0.270260183572876;
        case 22:
            return 0.264189241982359;
        case 23:
            return 0.2585097408088374;
        case 24:
            return 0.2531813568997612;
        case 25:
            return 0.2481693511764841;
        case 26:
            return 0.2434436124036166;
        case 27:
            return 0.2389778937254978;
        case 28:
            return 0.2347491976749156;
        case 29:
            return 0.2307372767004782;
        case 30:
            return 0.2269242244190854;
        case 31:
            return 0.2232941387424428;
        case 32:
            return 0.2198328424059889;
        case 33:
            return 0.216527649689546;
        case 34:
            return 0.2133671705705185;
        case 35:
            return 0.210341145412814;
        case 36:
            return 0.2074403047213377;
        case 37:
            return 0.204656249590948;
        case 38:
            return 0.2019813493339337;
        case 39:
            return 0.1994086534480554;
        case 40:
            return 0.1969318156005856;
        case 41:
            return 0.1945450278;
        case 42:
            return 0.1922429628;
        case 43:
            return 0.1900207248;
        case 44:
            return 0.1878738046;
        case 45:
            return 0.185798042;
        case 46:
            return 0.1837895915;
        case 47:
            return 0.1818448922;
        case 48:
            return 0.1799606416;
        case 49:
            return 0.1781337719;
        case 50:
            return 0.1763614288;
        case 51:
            return 0.1746409529;
        case 52:
            return 0.1729698629;
        case 53:
            return 0.1713458405;
        case 54:
            return 0.1697667173;
        case 55:
            return 0.1682304616;
        case 56:
            return 0.1667351687;
        case 57:
            return 0.1652790502;
        case 58:
            return 0.1638604244;
        case 59:
            return 0.1624777107;
        case 60:
            return 0.1611294174;
        case 61:
            return 0.1598141414;
        case 62:
            return 0.1585305558;
        case 63:
            return 0.1572774097;
        case 64:
            return 0.1560535159;
        default:
            std::cerr << "Ball <-> Cylinder mapping not yet implemented in dimension " << dim << std::endl;
            exit(1);
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
