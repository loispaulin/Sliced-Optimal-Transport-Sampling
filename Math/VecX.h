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

//VecXDynamic Declaration
class VecXDynamic {
private:
    int N;
    std::vector<double> coefs;
public:

    int dim() const;

    inline VecXDynamic(int n, const std::vector<double> &values);

    inline const double *data() const;

    inline double *data();

    explicit inline VecXDynamic(int n = 0);

    explicit inline operator double() const;

    inline VecXDynamic(const VecXDynamic &b);

    template<int n>
    inline VecXDynamic(const VecX<n> &v);

    inline VecXDynamic(double x, double y);

    inline VecXDynamic(double x, double y, double z);

    inline VecXDynamic(double x, double y, double z, double w);


    inline double norm() const;

    inline double norm2() const;

    inline void normalize();

    inline double &operator[](int i);

    inline const double &operator[](int i) const;

    inline VecXDynamic &operator+=(const VecXDynamic &b);

    inline VecXDynamic &operator-=(const VecXDynamic b);

    inline VecXDynamic &operator/=(double s);

    inline VecXDynamic &operator*=(double s);

    inline VecXDynamic &operator+=(double s);

    inline VecXDynamic &operator-=(double s);

    inline VecXDynamic operator-(const VecXDynamic &b) const;

    inline VecXDynamic operator+(const VecXDynamic &b) const;

    inline VecXDynamic operator/(const VecXDynamic &b) const;

    inline VecXDynamic operator-(double s) const;

    inline VecXDynamic operator+(double s) const;

    inline VecXDynamic operator*(double s) const;

    inline VecXDynamic operator/(double s) const;

    inline double operator*(const VecXDynamic &b) const;

    inline VecXDynamic operator-() const;

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

//VecXDynamic Definition
inline int VecXDynamic::dim() const {
    return N;
}

inline VecXDynamic::VecXDynamic(int n, const std::vector<double> &values) : N(n), coefs(values) {
    if (values.size() != size_t(N)) {
        std::cerr << "input vector must be of size " << N << " when creating VecXDynamic<" << N << std::endl;
        exit(1);
    }
}

inline const double *VecXDynamic::data() const {
    return coefs.data();
}

inline double *VecXDynamic::data() {
    return coefs.data();
}

inline VecXDynamic::VecXDynamic(int n) : N(n), coefs(n) {
}

inline VecXDynamic::operator double() const {
    if (N != 1) {
        std::cerr << "Error: Can't convert Vector to double if dimension is higher than 1." << std::endl;
        std::cerr << "\tCurrent dimension is " << N << std::endl;
        exit(1);
    }
    return coefs[0];
}

inline VecXDynamic::VecXDynamic(const VecXDynamic &b) : N(b.N), coefs(b.coefs) {}

template<int n>
inline VecXDynamic::VecXDynamic(const VecX<n> &v):N(n), coefs(n) {
    for (int i = 0; i < n; ++i) {
        coefs[i] = v[i];
    }
}

inline VecXDynamic::VecXDynamic(double x, double y) : N(2), coefs(2) {
    coefs[0] = x;
    coefs[1] = y;
}

inline VecXDynamic::VecXDynamic(double x, double y, double z) : N(3), coefs(3) {
    coefs[0] = x;
    coefs[1] = y;
    coefs[2] = z;
}

inline VecXDynamic::VecXDynamic(double x, double y, double z, double w) : N(4), coefs(4) {
    coefs[0] = x;
    coefs[1] = y;
    coefs[2] = z;
    coefs[3] = w;
}


inline double VecXDynamic::norm() const {
    return std::sqrt(norm2());
}

inline double VecXDynamic::norm2() const {
    double n = 0.;
    for (const double &v : coefs) {
        n += v * v;
    }
    return n;
}

inline void VecXDynamic::normalize() {
    double n = norm();
    for (double &v : coefs) {
        v /= n;
    }
}

inline double &VecXDynamic::operator[](int i) {
    return coefs[i];
}

inline const double &VecXDynamic::operator[](int i) const {
    return coefs[i];
}

inline VecXDynamic &VecXDynamic::operator+=(const VecXDynamic &b) {
    if (b.N != N) {
        std::cerr << "Operator +=  should be used between same sized VecXDynamic" << std::endl;
        exit(1);
    }
    for (size_t i = 0; i < coefs.size(); ++i) {
        coefs[i] += b.coefs[i];
    }
    return *this;
}

inline VecXDynamic &VecXDynamic::operator-=(const VecXDynamic b) {
    if (b.N != N) {
        std::cerr << "Operator -= should be used between same sized VecXDynamic" << std::endl;
        exit(1);
    }
    for (size_t i = 0; i < coefs.size(); ++i) {
        coefs[i] -= b.coefs[i];
    }
    return *this;

}

inline VecXDynamic &VecXDynamic::operator/=(double s) {
    for (size_t i = 0; i < coefs.size(); ++i) {
        coefs[i] /= s;
    }
    return *this;
}

inline VecXDynamic &VecXDynamic::operator*=(double s) {
    for (size_t i = 0; i < coefs.size(); ++i) {
        coefs[i] *= s;
    }
    return *this;
}

inline VecXDynamic &VecXDynamic::operator+=(double s) {
    for (size_t i = 0; i < coefs.size(); ++i) {
        coefs[i] += s;
    }
    return *this;
}

inline VecXDynamic &VecXDynamic::operator-=(double s) {
    for (size_t i = 0; i < coefs.size(); ++i) {
        coefs[i] -= s;
    }
    return *this;
}

inline VecXDynamic VecXDynamic::operator-(const VecXDynamic &b) const {
    if (b.N != N) {
        std::cerr << "Operator - should be used between same sized VecXDynamic" << std::endl;
        exit(1);
    }
    VecXDynamic a = *this;
    for (size_t i = 0; i < coefs.size(); ++i) {
        a.coefs[i] -= b.coefs[i];
    }
    return a;
}

inline VecXDynamic VecXDynamic::operator+(const VecXDynamic &b) const {
    if (b.N != N) {
        std::cerr << "Operator + should be used between same sized VecXDynamic" << std::endl;
        exit(1);
    }
    VecXDynamic a = *this;
    for (size_t i = 0; i < coefs.size(); ++i) {
        a.coefs[i] += b.coefs[i];
    }
    return a;
}

inline VecXDynamic VecXDynamic::operator/(const VecXDynamic &b) const {
    if (b.N != N) {
        std::cerr << "Operator / should be used between same sized VecXDynamic" << std::endl;
        exit(1);
    }
    VecXDynamic a = *this;
    for (size_t i = 0; i < coefs.size(); ++i) {
        a.coefs[i] /= b.coefs[i];
    }
    return a;
}

inline VecXDynamic VecXDynamic::operator-(double s) const {
    VecXDynamic a = *this;
    for (size_t i = 0; i < coefs.size(); ++i) {
        a.coefs[i] -= s;
    }
    return a;
}

inline VecXDynamic VecXDynamic::operator+(double s) const {
    VecXDynamic a = *this;
    for (size_t i = 0; i < coefs.size(); ++i) {
        a.coefs[i] += s;
    }
    return a;
}

inline VecXDynamic VecXDynamic::operator*(double s) const {
    VecXDynamic a = *this;
    for (size_t i = 0; i < coefs.size(); ++i) {
        a.coefs[i] *= s;
    }
    return a;
}

inline VecXDynamic VecXDynamic::operator/(double s) const {
    VecXDynamic a = *this;
    for (size_t i = 0; i < coefs.size(); ++i) {
        a.coefs[i] /= s;
    }
    return a;
}

inline double VecXDynamic::operator*(const VecXDynamic &b) const {
    if (b.N != N) {
        std::cerr << "Operator * should be used between same sized VecXDynamic" << std::endl;
        exit(1);
    }
    double res = 0;
    for (size_t i = 0; i < coefs.size(); ++i) {
        res += coefs[i] * b.coefs[i];
    }
    return res;
}

inline VecXDynamic VecXDynamic::operator-() const {
    VecXDynamic ret(*this);
    for (double &v : ret.coefs) {
        v = -v;
    }
    return ret;
}


inline VecXDynamic operator-(double s, const VecXDynamic &v) {
    VecXDynamic a(v.dim());
    for (int i = 0; i < a.dim(); ++i) {
        a[i] = s - v[i];
    }
    return a;
}

inline VecXDynamic operator+(double s, const VecXDynamic &v) {
    VecXDynamic a(v.dim());
    for (int i = 0; i < a.dim(); ++i) {
        a[i] = s + v[i];
    }
    return a;
}

inline VecXDynamic operator*(double s, const VecXDynamic &v) {
    VecXDynamic a(v.dim());
    for (int i = 0; i < a.dim(); ++i) {
        a[i] = s * v[i];
    }
    return a;
}

inline VecXDynamic operator/(double s, const VecXDynamic &v) {
    VecXDynamic a(v.dim());
    for (int i = 0; i < a.dim(); ++i) {
        a[i] = s / v[i];
    }
    return a;
}


inline bool operator<(const VecXDynamic &a, const VecXDynamic &b) {
    if (a.dim() != b.dim()) {
        std::cerr << "Operator < should be used between same sized VecXDynamic" << std::endl;
        exit(1);
    }
    for (int i = 0; i < a.dim(); ++i) {
        if (a[i] < b[i])
            return true;
        else if (a[i] > b[i])
            return false;
    }
    return false;
}


inline double scalar(const VecXDynamic &a, const VecXDynamic &b) {
    return a * b;
}


inline std::ostream &operator<<(std::ostream &out, const VecXDynamic &v) {
    for (int i = 0; i < v.dim() - 1; ++i) {
        out << v[i] << " ";
    }
    out << v[v.dim() - 1];
    return out;
}


#endif //SLICEDOPTIM_VECN_H
