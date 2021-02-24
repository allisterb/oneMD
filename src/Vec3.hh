
#ifndef VECTOR_H
#define VECTOR_H

#define CHUNKSIZE 15
#include "xdrfile.h"
#include <array>
using namespace std;

/** X coordinate */
const int X = 0;
/** Y coordinate */
const int Y = 1;
/** Z coordinate */
const int Z = 2;

class Vec3 {

private:
    array <double,3> r;
public:
    Vec3();
    Vec3(double x, double y, double z);
    void set(double x, double y, double z);
    Vec3 operator-(Vec3 rhs);
    void operator -= (Vec3 rhs);
    Vec3 operator + (Vec3 rhs);
    void operator += (Vec3 rhs);
    Vec3 operator / (double rhs);
    void operator /= (double rhs);
    void operator *= (double rhs);
    void operator = (double rhs);
    friend Vec3 operator * (Vec3 lhs, double rhs);
    friend Vec3 operator*(double lhs, Vec3 rhs);
    double& operator[](int i);
    const double& operator[](int i) const;
};

using Vec3 = Vec3;
#endif
