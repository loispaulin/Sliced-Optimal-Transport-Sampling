#include <sstream>
#include "iopointset.h"


bool read_points_from_file(std::istream& in, std::vector<double>& points){
    points.clear();
    std::string line;
    while(std::getline(in, line)){
        int c = line.find_first_not_of(" \t");
        if (line[c] != '#'){
            std::istringstream lineIn(line);
            double d;
            while (lineIn >> d){
                points.push_back(d);
            }
        } else {
            return true;
        }
    }
    return points.size() != 0;
}