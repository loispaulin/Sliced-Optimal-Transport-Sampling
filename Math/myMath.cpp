//
//

#include <cmath>
#include <vector>
#include <iostream>
#include "myMath.h"


long int binomialCoef(int n, int k){
    if (k > n / 2){
        k = n - k;
    }
    long int res = 1;
    for (int i = 1; i <= k; ++i){
        res *= n - k + i;
        res /= i;
    }
    return res;
}

/**
 * Return k among n using the value of (k-1) among n
 * @param n
 * @param k
 * @param prev (k-1) among n
 * @return k among n
 */
long int binomialCoefNext(int n, int k, long int prev){
    return prev * (n - k + 1) / k;
}
/**
 * Return k among n using the value of (k-1) among n
 * @param n
 * @param k
 * @param prev (k-1) among n
 * @return k among n
 */
double binomialCoefNext(int n, int k, double prev){
    return prev * (n - k + 1) / k;
}

double myfmod(double v, double m){
    return std::fmod(std::fmod(v, m) + m, m);
}

double nballVolume(double r, int dim){
    double res = std::pow(r, dim);
    res *= std::pow(M_PI, dim / 2.);
    if (dim % 2){
        for (int i = 1; i <= (dim) / 2 + 1; ++i){
            res /= (i - 0.5);
        }
        res /= std::sqrt(M_PI);
    } else {
        for (int i = 1; i <= dim / 2; ++i){
            res /= i;
        }
    }
    return res;
}

double nsphereArea(double r, int dim){
    double res = 2 * std::pow(M_PI, (dim) * 0.5);
    res /= tgamma((dim) * 0.5);
    res *= std::pow(r, dim-1);
    return res;
}


int getIntLog2(int v){
    int res = 0;
    while (v){
        res += 1;
        v /= 2;
    }
    return res - 1;
}


void neumaierStep(double& res, double& val, double& comp){

    double t = res + val;
    if (std::abs(res) > std::abs(val)){
        comp += (res - t) + val;
    } else {
        comp += (val - t) + res;
    }
    res = t;

}

double integrateSinPowerOdd(double a, double b, int power){

    int n = power / 2;

    double res = 0.;
    double comp = 0.;
    int currentFact = n;
    double factor = 1.;
    static const double maxi = std::pow(10., 5.);
    static const double mini = std::pow(10., -4.);

    for (int i = 0; i <= n; ++i){
        double insidefactor = factor / (2. * i + 1.);
        double val = insidefactor * (std::pow(std::cos(b), 2. * i + 1.) - std::pow(std::cos(a), 2. * i + 1.));
        neumaierStep(res, val, comp);
        factor = -binomialCoefNext(n, i+1, factor);
        while(currentFact > 1 && (std::abs(factor) > maxi || std::abs(res) > maxi)){
            factor /= currentFact;
            res /= currentFact;
            comp /= currentFact;
            currentFact -= 1;
        }
        while(currentFact < n && (std::abs(factor) < mini || std::abs(res) < mini)){
            factor *= currentFact + 1;
            res *= currentFact + 1;
            comp *= currentFact + 1;
            currentFact += 1;
        }
    }
    for (int i = currentFact; i < n; ++i){
        res *= currentFact + 1;
        comp *= currentFact + 1;
    }

    return -(res + comp);
}

double integrateCosPowerOdd(double a, double b, int power){

    int n = power / 2;

    double res = 0.;
    double comp = 0.;
    int currentFact = n;
    double factor = 1.;
    static const double maxi = std::pow(10., 6.);
    static const double mini = std::pow(10., -3.);

    for (int i = 0; i <= n; ++i){
        double insidefactor = factor / (2. * i + 1.);
        double val = insidefactor * (std::pow(std::sin(b), 2 * i + 1) - std::pow(std::sin(a), 2 * i + 1));
        neumaierStep(res, val, comp);

        factor = -binomialCoefNext(n, i+1, factor);

        while(currentFact > 1 && (std::abs(factor) > maxi || std::abs(res) > maxi)){
            factor /= currentFact;
            res /= currentFact;
            comp /= currentFact;
            currentFact -= 1;
        }
        while(currentFact < n && (std::abs(factor) < mini || std::abs(res) < mini)){
            factor *= currentFact + 1;
            res *= currentFact + 1;
            comp *= currentFact + 1;
            currentFact += 1;
        }
    }
    for (int i = currentFact; i < n; ++i){
        res *= currentFact + 1;
        comp *= currentFact + 1;
    }

    return res + comp;
}

/**
 * Computes \$f \int_a^b sin(x)^{p} dx \$f numerical precision make it so that it doesn't work for p > 120
 * @param a
 * @param b
 * @param power p
 * @return
 */
double integrateSinPower(double a, double b, int power){

    //If we are in an imprecise case use numerical integration instead of closed form formula
    if (power > 32 && std::abs(std::fmod(b, M_PI) - M_PI_2) > 0.05){
        std::function<double(double)> f = [power](double x){return std::pow(std::sin(x), power); };
        return simpsonIntegration(f, a, b, 1000);
    }

    //Cos^0 = 1 and \int_a^b 1 dx = b - a
    if (power == 0){
        //std::cout << "{\"NIntegrate[Sin[x]^" << power << ", {x, " << a << ", " << b << "} ]\", NIntegrate[Sin[x]^" << power << ", {x, " << a << ", " << b << "} ] - " << b - a << " < 0.000001}, " << std::endl;
        return b - a;
    }
    //Cos^{2*n+1} has a closed form
    if (power % 2){
        double res = integrateSinPowerOdd(a, b, power);
        //std::cout << "{\"NIntegrate[Sin[x]^" << power << ", {x, " << a << ", " << b << "} ]\", NIntegrate[Sin[x]^" << power << ", {x, " << a << ", " << b << "} ] - " << res << " < 0.000001}, " << std::endl;
        return res;
    }

    //Else dynamic computation
    int tabSize = 1 + (power / 2);
    std::vector<double> T(2 * tabSize);
T.data();
    int nbSteps = getIntLog2(power);

    for (int i = nbSteps; i > 0; --i){

        for (int j = 0; j <= power / std::pow(2., i); ++j){
            //Computing Integral[ Cos[x]^j, {x, 2^i * a, 2^i * b} ]
            if (j == 0){

                T[(i%2) * tabSize + j] = std::pow(2., i) * b - std::pow(2., i) * a;

            } else if (j % 2){

                T[(i%2) * tabSize + j] = integrateCosPowerOdd(std::pow(2., i) * a, std::pow(2., i) * b, j);

            } else {

                double res = 0.;
                double comp = 0.;
                int n = j / 2;
                int currentPow = 0;
                double factor = 1.;
                static const double maxi = std::pow(10., 6.);


                for (int k = 0; k <= n; ++k){
                    double val = factor * T[((i+1)%2) * tabSize + k];
                    neumaierStep(res, val, comp);
                    factor = binomialCoefNext(n, k+1, factor);
                    if (currentPow < n + 1 && (res > maxi || factor > maxi)){
                        res /= 2.;
                        factor /= 2.;
                        comp /= 2.;
                        currentPow += 1;
                    }
                }
                double remainingPow = std::pow(2., n + 1 - currentPow);
                res /= remainingPow;
                comp /= remainingPow;

                T[(i%2) * tabSize + j] = res + comp;

            }

        }

    }

    double res = 0.;
    double comp = 0.;
    int n = power / 2;
    int currentPow = 0;
    double factor = 1.;
    static const double maxi = std::pow(10., 6.);

    for (int k = 0; k <= n; ++k){
        double val = factor * T[tabSize + k];
        neumaierStep(res, val, comp);
        factor = -binomialCoefNext(n, k+1, factor);
        if (currentPow < n + 1 && (std::abs(res) > maxi || std::abs(factor) > maxi)){
            res /= 2.;
            comp /= 2.;
            factor /= 2.;
            currentPow += 1;
        }
    }
    double remainingPow = std::pow(2., n + 1 - currentPow);
    res /= remainingPow;
    comp /= remainingPow;

    //std::cout << "{\"NIntegrate[Sin[x]^" << power << ", {x, " << a << ", " << b << "} ]\", NIntegrate[Sin[x]^" << power << ", {x, " << a << ", " << b << "} ] - " << res + comp << " < 0.000001}, " << std::endl;

    return res + comp;

}


/**
 * Computes \$f \int_a^b cos(x)^{p} dx \$f numerical precision make it so that it doesn't work for p > 120
 * @param a
 * @param b
 * @param power p
 * @return
 */
double integrateCosPower(double a, double b, int power){


    //Cos^0 = 1 and \int_a^b 1 dx = b - a
    if (power == 0){
        //std::cout << "{\"NIntegrate[Cos[x]^" << power << ", {x, " << a << ", " << b << "} ]\", NIntegrate[Cos[x]^" << power << ", {x, " << a << ", " << b << "} ] - " << b - a << " < 0.000001}, " << std::endl;
        return b - a;
    }
    //Cos^{2*n+1} has a closed form
    if (power % 2){
        double res = integrateCosPowerOdd(a, b, power);
        //std::cout << "{\"NIntegrate[Cos[x]^" << power << ", {x, " << a << ", " << b << "} ]\", NIntegrate[Cos[x]^" << power << ", {x, " << a << ", " << b << "} ] - " << res << " < 0.000001}, " << std::endl;
        return res;
    }

    //Else dynamic computation
    int tabSize = 1 + (power / 2);
    std::vector<double> T(2 * tabSize);

    int nbSteps = getIntLog2(power);

    for (int i = nbSteps; i > 0; --i){

        for (int j = 0; j <= power / std::pow(2., i); ++j){

            if (j == 0){

                T[(i%2) * tabSize + j] = std::pow(2., i) * b - std::pow(2., i) * a;

            } else if (j % 2){

                T[(i%2) * tabSize + j] = integrateCosPowerOdd(std::pow(2., i) * a, std::pow(2., i) * b, j);

            } else {

                double res = 0.;
                double comp = 0.;
                int n = j / 2;
                int currentPow = 0;
                double factor = 1.;
                static const double maxi = std::pow(10., 9.);

                for (int k = 0; k <= n; ++k){
                    double val = factor * T[((i+1)%2) * tabSize + k];
                    neumaierStep(res, val, comp);
                    factor = binomialCoefNext(n, k+1, factor);
                    if (currentPow < n + 1 && (std::abs(res) > maxi || std::abs(factor) > maxi)){
                        res /= 2.;
                        comp /= 2.;
                        factor /= 2.;
                        currentPow += 1;
                    }
                }
                double remainingPow = std::pow(2., n + 1 - currentPow);
                res /= remainingPow;
                comp /= remainingPow;

                T[(i%2) * tabSize + j] = res + comp;

            }

        }

    }

    double res = 0.;
    double comp = 0.;
    int n = power / 2;
    int currentPow = 0;
    double factor = 1.;
    static const double maxi = std::pow(10., 9.);

    for (int k = 0; k <= n; ++k){
        double val = factor * T[tabSize + k];
        neumaierStep(res, val, comp);
        factor = binomialCoefNext(n, k+1, factor);
        if (currentPow < n + 1 && (res > maxi || factor > maxi)){
            res /= 2.;
            factor /= 2.;
            currentPow += 1;
        }
    }
    double remainingPow = std::pow(2., n + 1 - currentPow);
    res /= remainingPow;
    comp /= remainingPow;

    //std::cout << "{\"NIntegrate[Cos[x]^" << power << ", {x, " << a << ", " << b << "} ]\", NIntegrate[Cos[x]^" << power << ", {x, " << a << ", " << b << "} ] - " << res + comp << " < 0.000001}, " << std::endl;

    return res + comp;

}

/**
 * Evaluate the integral of \param f between \param a and \param b using Simpson's rule with \param n + 1 points
 * @param f Function to integrate
 * @param a
 * @param b
 * @param n
 * @return
 */
double simpsonIntegration(std::function<double(double)> f, double a, double b, int n){
    if (n % 2){
        std::cerr << "Simpson Integration error: n must be even" << std::endl;
        exit(1);
    }
    double dx = (b - a) / n;
    double res = 0.;

    for (int i = 0; i < n / 2; ++i){
        res += 2 * f(a + 2 * i * dx) + 4 * f(a + (2 * i + 1) * dx);
    }

    res += f(b);
    //we added an extra f(a) in the loop (2 * f(a + 0 * dx)) so we remove it
    res -= f(a);

    res *= dx / 3.;

    return res;
}




















