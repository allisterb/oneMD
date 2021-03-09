/* Based on https://github.com/wesbarnett/lennardjones" by James W. Barnett
 * libgmxcpp
 * Copyright (C) 2015 James W. Barnett <jbarnet4@tulane.edu>
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * The full license is located in a text file titled "LICENSE" in the root
 * directory of the source.
 *
 */

#pragma once

#include <array>
#include "xdrfile.h"

#include "common.hh"

using namespace std;

/** X coordinate */
constexpr int X = 0;
/** Y coordinate */
constexpr int Y = 1;
/** Z coordinate */
constexpr int Z = 2;

#ifndef USE_ONEAPI
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
#else
using Vec3 = sycl::double3;
#endif