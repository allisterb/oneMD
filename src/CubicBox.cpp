/* Based on libgmxcpp by James W. Barnett: https://github.com/wesbarnett/libgmxcpp
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

#include "CubicBox.hh"

#ifndef USE_ONEAPI
#include <cmath>

CubicBox::CubicBox() {}

CubicBox::CubicBox(float x, float y, float z)
{
    this->box[X] = x;
    this->box[Y] = y;
    this->box[Z] = z;
}
float& CubicBox::operator[](int i)
{
    return box[i];
}

const float& CubicBox::operator[](int i) const
{
    return box[i];
}

Vec3 pbc(Vec3 a, CubicBox box)
{
    a[Z] -= box[Z] * nearbyint(a[Z] / box[Z]);
    a[Y] -= box[Y] * nearbyint(a[Y] / box[Y]);;
    a[X] -= box[X] * nearbyint(a[X] / box[X]);;
    return a;
}
double distance(Vec3 a, Vec3 b, CubicBox box)
{
    return sqrt(distance2(a, b, box));
}

double distance2(Vec3 a, Vec3 b, CubicBox box)
{
    Vec3 c = a - b;
    c = pbc(c, box);
    return dot(c, c);
}
double magnitude(Vec3 x)
{
    return sqrt(dot(x, x));
}

double volume(CubicBox box)
{
    return box[X] * box[Y] * box[Z];
}

double dot(Vec3 a, Vec3 b)
{
    return a[X] * b[X] + a[Y] * b[Y] + a[Z] * b[Z];
}
#else
Vec3 pbc(Vec3 a, CubicBox box)
{
    return a - box * sycl::rint(a / box);
}

double distance(Vec3 a, Vec3 b, CubicBox box)
{
    return sycl::sqrt(distance2(a, b, box));
}

double distance2(Vec3 a, Vec3 b, CubicBox box)
{
    Vec3 c = a - b;
    c = pbc(c, box);
    return dot(c, c);
}

double magnitude(Vec3 x)
{
    return sycl::sqrt(dot(x, x));
}

double volume(CubicBox box)
{
    return box[X] * box[Y] * box[Z];
}

double dot(Vec3 a, Vec3 b)
{
    return sycl::dot(a, b);
}
#endif