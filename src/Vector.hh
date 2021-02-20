
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

class Vector {

private:
    array <double,3> r;
public:
    Vector();
    Vector(double x, double y, double z);
    void set(double x, double y, double z);
    Vector operator-(Vector rhs);
    void operator -= (Vector rhs);
    Vector operator + (Vector rhs);
    void operator += (Vector rhs);
    Vector operator / (double rhs);
    void operator /= (double rhs);
    void operator *= (double rhs);
    void operator = (double rhs);
    friend Vector operator * (Vector lhs, double rhs);
    friend Vector operator*(double lhs, Vector rhs);
    double& operator[](int i);
    const double& operator[](int i) const;
};

using Vec3 = Vector;
#endif
