//
//

#include <cmath>
#include <sstream>
#include <fstream>
#include <cstring>
#include "my_utility.h"

using namespace std;

double clamp(double v, double min, double max){
    return v > min ? (v < max ? v : max) : min;
}

double inverseFunction(std::function<double(double)>& f, std::function<double(double)>& df, double y){
    static const double eps = std::pow(10, -14);

    double x = 0.;
    double act = f(x);

    int nbsteps = 0;
    while (std::abs(y - act) > eps) {
        double g = df(x);
        double s =  (y - act) / g;
        x += s;
        if (abs(s) < eps){
            break;
        }
        act = f(x);
        ++nbsteps;
        if (nbsteps > 100) {
            break;
        }
    }

    return x;
}
