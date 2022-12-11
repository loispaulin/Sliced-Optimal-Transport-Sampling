//
// Created by lpaulin on 03/12/19.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstring>
#include <sstream>
#include <random>
#include <omp.h>
#include "../Math/myMath.h"
#include "../Math/VecX.h"
#include "../Tools/iopointset.h"
#include "../Tools/mapping.h"
#include "../Transport/slicedOptimalTransportNBall.h"

#define DIM 2

using namespace std;


template <typename T>
inline std::istream& operator>>(std::istream& in, std::vector<T>& v){

    v.clear();
    char c;
    in >> c;
    if (c != '{'){
        std::cerr << "Error: Incorrect input format for operator >> on a vector '{' expected" << std::endl;
        exit(1);
    }
    while (c != '}'){
        T val;
        in >> val;
        v.push_back(std::move(val));
        in >> c;
    }

    return in;

}

void usage(const char **argv){
    cerr << argv[0] << "[-i <inputFileName>]  [-o <OutputFileName>] "
                       "[-n <nbPoints>] [-m <nbRealisations>] [-p <nbIteration>] "
                       "[--step <nbDirectionPerStep>] [-s <seed>] [-d <dimension>] [-c] "
                       "[--proj <projections> --sizes <nbSlicesPerDim>]" << endl;
}

void handleParameters(int argc,
                      const char** argv,
                      string& filename,
                      string& outPrefix,
                      int& nbIter,
                      int& m,
                      int& p,
                      int& seed,
                      int& nbPointsets,
                      int& dim,
                      bool& cubify,
                      vector<vector<bool>>& projections,
                      vector<size_t>& projSizes,
                      vector<double>& weights,
                      bool& silent,
                      bool& initFromCube){
    int i = 1;
    while (i < argc){
        if (!strncmp(argv[i], "-j", 2)){
            omp_set_num_threads(atoi(argv[i+1]));
            ++i;
        } else if (!strncmp(argv[i], "-i", 2)){
            filename = argv[i+1];
            ++i;
        } else if (!strncmp(argv[i], "-o", 2)) {
            outPrefix = (argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "-n", 2)) {
            p = atoi(argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "--step", 6)) {
            m = atoi(argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "-m", 2)) {
            nbPointsets = atoi(argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "-p", 2)) {
            nbIter = atoi(argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "-s", 2)) {
            seed = atoi(argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "-d", 2)) {
            dim = atoi(argv[i+1]);
            ++i;
        } else if (!strncmp(argv[i], "-c", 2)) {
            cubify = true;
        } else if (!strncmp(argv[i], "--fromCube", 10)) {
            initFromCube = true;
        } else if (!strncmp(argv[i], "--proj", 6)) {
            std::istringstream in(argv[i+1]);
            vector<vector<size_t>> indices;
            in >> indices;
            projections = vector<vector<bool>>(indices.size(), vector<bool>(dim, false));
            for(size_t ind = 0; ind < indices.size(); ++ind){
                for (const size_t& v : indices[ind]){
                    projections[ind][v] = true;
                }
            }
            ++i;
        } else if (!strncmp(argv[i], "-w", 2)) {
            std::istringstream in(argv[i+1]);
            in >> weights;
            double sum = 0.;
            for(double d : weights){
                sum += d;
            }
            for(double& d : weights){
                d /= sum;
            }
            ++i;
        } else if (!strncmp(argv[i], "--sizes", 7)) {
            std::istringstream in(argv[i+1]);
            in >> projSizes;
            for(size_t ind = 1; ind < projSizes.size(); ++ind){
                projSizes[ind] += projSizes[ind-1];
            }
            ++i;
        } else if (!strncmp(argv[i], "-h", 2) || !strncmp(argv[i], "--help", 6)) {
            cerr << "Help: " << endl;
            cerr << "Option list:" << endl;
            cerr << "\t-i <inputFileName> (optional): Specifies an input file in which to read points."
                 << "As many will be read as the -p parameter specifies" << endl;
            cerr << "\t-o <OutputFileName> (optional): Specifies an output file in which points will be written."
                 << "If unset standard output will be used" << endl;
            cerr << "\t-n <nbPoints> (default 1024): Specifies the number of points to generate" << endl;
            cerr << "\t-m <nbRealisations> (default 1): Specifies the number of generated pointsets" << endl;
            cerr << "\t-p <nbIteration> (default 1024): Specifies the number of batches in the optimization process"
                    << endl;
            cerr << "\t--step <nbDirectionPerStep> (default 32): Specifies the number of slices per batch in "
                    "the optimization process" << endl;
            cerr << "\t-s <seed> (default 133742): Specifies the random seed" << endl;
            cerr << "\t-d <dimension> (default 2): Specifies the dimension of the points to generate" << endl;
            cerr << "\t-c (optional): If unset points will be given in the unit ball. Else they will be in "
                    "the unit cube [0,1)^d" << endl;
            cerr << "\t--proj (optional): C style definition of an array of arrays of integers representing which dimensions "
                    "are active in each projection (sizes must be set too). WARNING: '-d' must be used before." << endl;
            cerr << "\t--sizes (optional): C style definition of a array of doubles representing the number of slices "
                    "to use for each projection" << endl;
            cerr << "\t--silent (optional): Cancels all outputs other than the points and errors" << endl;
            cerr << "\t--fromCube (optional): Maps input initialisation pointset from cube to ball" << endl;
            cerr << "\t" << endl;
            usage(argv);
            exit(2);
        } else if (!strncmp(argv[i], "--silent", 8)) {
            silent = true;
        } else{
            cerr << "Unknown parameters! Use -h to get the option list." << std::endl;
            exit(3);
        }
        ++i;
    }
}

template <class VECTYPE>
int main_template(int argc, const char **argv) {

    int nbIter = 512;
    int dim = 4;
    int p = 4;
    //Default parameters value
    string filename;    //int nbIter = 2;
    string outPrefix;
    int seed = 133742;
    int nbPointsets = 1024;
    bool cubify = false;
    bool silent = false;
    bool initFromCube = false;
    vector<vector<size_t>> bla = { {0,1,2,3,4,5},{0,1}, {2,3}, {4,5} };
    vector<vector<bool>> projections = vector<vector<bool>>(bla.size(), vector<bool>(dim, false));// = {{1,1,1,1,1,1}};
    vector<double> projWeight(projections.size(), 1. / double(projections.size()));
    for(size_t ind = 0; ind < bla.size(); ++ind){
        for (const size_t& v : bla[ind]){
            projections[ind][v] = true;
        }
    }

    int m = 80;
    vector<size_t> projSizes = {32,16,16,16};
    for(size_t ind = 1; ind < projSizes.size(); ++ind){
        projSizes[ind] += projSizes[ind-1];
    }
    //Read parameters from call
    handleParameters(argc, argv, filename, outPrefix, nbIter, m, p, seed, nbPointsets, dim, cubify, projections,
                     projSizes, projWeight, silent, initFromCube);

    if (projections.size() != projSizes.size()){
        cerr << "Error: Projections array and projections sizes array are of different size." << endl;
        exit(1);
    }

    if (projSizes.back() != size_t(m)){
        cerr << "Error: Wrong number of slices per projection. They don't add up to " << m << endl;
        exit(1);
    }

    int origP = p;

    //If file name ends in .bin then the output will be written in binary
    ostream* out = &cout;
    if (outPrefix != "") {
        out = new ofstream(outPrefix);
        if (out->fail()) {
            cerr << "Can't open output file \"" + outPrefix + "\"" << endl;
            exit(3);
        }
        cerr << "Output file: " << outPrefix << endl;
    }

    if (!silent) {
        cerr << "Generating " << nbPointsets << " sets of " << origP << " points in " << dim << "D using " << nbIter << " batches of " << m
             << " slices" << endl;
        cerr << "Map points to cube: " << (cubify ? "True" : "False") << endl;
    }

    for (int indPointset = 0; indPointset < nbPointsets; ++indPointset) {
        mt19937 generator(seed + indPointset);
        vector<VECTYPE> points(p, VECTYPE(dim));

        if (filename == "") {
            //Init from whitenoise in ball
            uniform_real_distribution<> unif(-1.0, 1.0);
            normal_distribution<> norm(0., 1.0);

            for (size_t i = 0; i < points.size(); ++i) {
                for (size_t j = 0; j < size_t(dim); ++j){
                    points[i][j] = unif(generator);
                }
            }

        } else {
            //Init from points in file
            std::ifstream in(filename);
            if (in.fail()) {
                cerr << "Error opening input file \"" + filename + "\"" << endl;
            }
            vector<double> linearpoints;
            read_points_from_file(in, linearpoints);

            points.resize(linearpoints.size() / dim);
            //Store points in pair vector
            for (size_t i = 0; i < points.size(); ++i) {
                for (int j = 0; j < dim; ++j) {
                    points[i][j] = linearpoints[dim * i + j];
                }
                if (initFromCube){
                    points[i] = NCube2NBall(2. * points[i] - 1.);
                }
            }
        }

        vector<VECTYPE> result;

        newSlicedOptimalTransport(points, result, nbIter, m, generator, projections, projSizes, projWeight);

        if (indPointset != 0)
            *out << "#" << endl;

        for (VECTYPE& v : result){
            v = 0.5 * v + 0.5;
        }
        savePointsetND(*out, result);

    }

    return 0;
}

#include "../Tools/dimensionsInstantiation.hpp"
