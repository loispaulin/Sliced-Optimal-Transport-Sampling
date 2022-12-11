//
//

#ifndef SLICEDOPTIM_IOPOINTSET_H
#define SLICEDOPTIM_IOPOINTSET_H

#include <iostream>
#include <fstream>
#include <vector>
#include "../Math/VecX.h"

bool read_points_from_file(std::istream& in, std::vector<double>& points);

template <class VECTYPE>
inline void savePointsetND(std::ostream& out, const std::vector<VECTYPE>& points){
    out.precision(11);
    for (auto &v : points) {
        out << v << std::endl;
    }
}

#endif //SLICEDOPTIM_IOPOINTSET_H
