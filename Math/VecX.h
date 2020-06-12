//
//

#ifndef SLICEDOPTIM_VECN_H
#define SLICEDOPTIM_VECN_H

#include <vector>
#include <iostream>
#include <cmath>

template<int N>
class VecX;

//VecX Declaration
template<int N>
class VecX {
private:
    double coefs[N];
public:

    inline VecX();

    explicit inline VecX(int dim);

    inline VecX(double x, double y);

    inline VecX(const std::vector<double> &values);

    const double *data() const;

    double *data();

    inline int dim() const;

    inline double norm() const;

    inline double norm2() const;

    inline void normalize();

    inline double &operator[](int i);

    inline const double &operator[](int i) const;

    inline VecX &operator+=(const VecX<N> &b);

    inline VecX &operator-=(const VecX b);

    inline VecX &operator/=(double s);

    inline VecX &operator*=(double s);

    inline VecX &operator+=(double s);

    inline VecX &operator-=(double s);

    inline VecX operator-(const VecX &b) const;

    inline VecX operator+(const VecX &b) const;

    inline VecX operator/(const VecX &b) const;

    inline VecX operator-(double s) const;

    inline VecX operator+(double s) const;

    inline VecX operator*(double s) const;

    inline VecX operator/(double s) const;

    inline double operator*(const VecX &b) const;

};

//VecX Definition

template<int N>
inline VecX<N>::VecX(const std::vector<double> &values) {
    if (values.size() != N) {
        std::cerr << "input vector must be of size " << N << " when creating VecX<" << N << std::endl;
    }
    memcpy(coefs, &values[0], N * sizeof(double));
}

template<int N>
inline VecX<N>::VecX() {}

template<int N>
inline VecX<N>::VecX(int dim) {
    if (dim != N) {
        std::cerr << "input vector must be of size " << N << " when creating VecX<" << N << "> but received " << dim
                  << std::endl;
    }
}

template<int N>
inline VecX<N>::VecX(double x, double y) {
    if (N < 2) {
        std::cerr << "Creating a Vec1 from 2 values" << std::endl;
        exit(1);
    }
    coefs[0] = x;
    coefs[1] = y;
}

template<int N>
inline const double *VecX<N>::data() const {
    return coefs;
}

template<int N>
inline double *VecX<N>::data() {
    return coefs;
}

template<int N>
inline int VecX<N>::dim() const { return (int) N; }

template<int N>
inline double VecX<N>::norm() const {
    return std::sqrt(norm2());
}

template<int N>
inline double VecX<N>::norm2() const {
    double n = 0.;
    for (int i = 0; i < N; i++) {
        double v = coefs[i];
        n += v * v;
    }
    return n;
}

template<int N>
inline void VecX<N>::normalize() {
    double n = norm();
    for (int i = 0; i < N; i++) {
        coefs[i] /= n;
    }
}

template<int N>
inline double &VecX<N>::operator[](int i) {
    return coefs[i];
}

template<int N>
inline const double &VecX<N>::operator[](int i) const {
    return coefs[i];
}

template<int N>
inline VecX<N> &VecX<N>::operator+=(const VecX<N> &b) {
    for (int i = 0; i < N; ++i) {
        coefs[i] += b.coefs[i];
    }
    return *this;
}

template<int N>
inline VecX<N> &VecX<N>::operator-=(const VecX<N> b) {
    for (int i = 0; i < N; ++i) {
        coefs[i] -= b.coefs[i];
    }
    return *this;

}

template<int N>
inline VecX<N> &VecX<N>::operator/=(double s) {
    for (int i = 0; i < N; ++i) {
        coefs[i] /= s;
    }
    return *this;
}

template<int N>
inline VecX<N> &VecX<N>::operator*=(double s) {
    for (int i = 0; i < N; ++i) {
        coefs[i] *= s;
    }
    return *this;
}

template<int N>
inline VecX<N> &VecX<N>::operator+=(double s) {
    for (int i = 0; i < N; ++i) {
        coefs[i] += s;
    }
    return *this;
}

template<int N>
inline VecX<N> &VecX<N>::operator-=(double s) {
    for (int i = 0; i < N; ++i) {
        coefs[i] -= s;
    }
    return *this;
}

template<int N>
inline VecX<N> VecX<N>::operator-(const VecX<N> &b) const {
    VecX<N>
            a = *this;
    for (int i = 0; i < N; ++i) {
        a.coefs[i] -= b.coefs[i];
    }
    return a;
}

template<int N>
inline VecX<N> VecX<N>::operator+(const VecX<N> &b) const {
    VecX<N>
            a = *this;
    for (int i = 0; i < N; ++i) {
        a.coefs[i] += b.coefs[i];
    }
    return a;
}

template<int N>
inline VecX<N> VecX<N>::operator/(const VecX<N> &b) const {
    VecX<N>
            a = *this;
    for (int i = 0; i < N; ++i) {
        a.coefs[i] /= b.coefs[i];
    }
    return a;
}

template<int N>
inline VecX<N> VecX<N>::operator-(double s) const {
    VecX<N>
            a = *this;
    for (int i = 0; i < N; ++i) {
        a.coefs[i] -= s;
    }
    return a;
}

template<int N>
inline VecX<N> VecX<N>::operator+(double s) const {
    VecX<N>
            a = *this;
    for (int i = 0; i < N; ++i) {
        a.coefs[i] += s;
    }
    return a;
}

template<int N>
inline VecX<N> VecX<N>::operator*(double s) const {
    VecX<N>
            a = *this;
    for (int i = 0; i < N; ++i) {
        a.coefs[i] *= s;
    }
    return a;
}

template<int N>
inline VecX<N> VecX<N>::operator/(double s) const {
    VecX<N>
            a = *this;
    for (int i = 0; i < N; ++i) {
        a.coefs[i] /= s;
    }
    return a;
}

template<int N>
inline double VecX<N>::operator*(const VecX<N> &b) const {
    double res = 0;
    for (int i = 0; i < N; ++i) {
        res += coefs[i] * b.coefs[i];
    }
    return res;
}

template<int N>
inline VecX<N> operator-(const VecX<N> &v) {
    VecX<N> v1;
    for (int i = 0; i < N; ++i) {
        v1[i] = -v[i];
    }
    return v1;
}

template<int N>
inline VecX<N> operator-(double s, const VecX<N> &v) {
    VecX<N> a;
    for (int i = 0; i < N; ++i) {
        a[i] = s - v[i];
    }
    return a;
}

template<int N>
inline VecX<N> operator+(double s, const VecX<N> &v) {
    VecX<N> a;
    for (int i = 0; i < N; ++i) {
        a[i] = s + v[i];
    }
    return a;
}

template<int N>
inline VecX<N> operator*(double s, const VecX<N> &v) {
    VecX<N> a;
    for (int i = 0; i < N; ++i) {
        a[i] = s * v[i];
    }
    return a;
}

template<int N>
inline VecX<N> operator/(double s, const VecX<N> &v) {
    VecX<N> a;
    for (int i = 0; i < N; ++i) {
        a[i] = s / v[i];
    }
    return a;
}

template<int N>
inline VecX<N> operator/(const VecX<N> &v, double s) {
    VecX<N> a;
    for (int i = 0; i < N; ++i) {
        a[i] = v[i] / s;
    }
    return a;
}

template<int N>
inline bool operator<(const VecX<N> &a, const VecX<N> &b) {
    for (int i = 0; i < N; ++i) {
        if (a[i] < b[i])
            return true;
        else if (a[i] > b[i])
            return false;
    }
}

template<int N>
inline double scalar(const VecX<N> &a, const VecX<N> &b) {
    return a * b;
}

template<int N>
inline std::ostream &operator<<(std::ostream &out, const VecX<N> &v) {
    for (int i = 0; i < N - 1; ++i) {
        out << v[i] << " ";
    }
    out << v[N - 1];
    return out;
}

template<int N>
inline VecX<N> mod(const VecX<N> &v, double m) {
    VecX<N> res;
    for (int i = 0; i < N; ++i) {
        res[i] = std::fmod(v[i], m);
    }
    return res;
}

#endif //SLICEDOPTIM_VECN_H
