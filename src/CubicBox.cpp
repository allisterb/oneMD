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
    a[Z] -= box[Z] * __nearbyint(a[Z] / box[Z]);
    a[Y] -= box[Y] * __nearbyint(a[Y] / box[Y]);;
    a[X] -= box[X] * __nearbyint(a[X] / box[X]);;
    return a;
}

double distance(Vec3 a, Vec3 b, CubicBox box)
{
    return __sqrt(distance2(a, b, box));
}

double distance2(Vec3 a, Vec3 b, CubicBox box)
{
    Vec3 c = a - b;
    c = pbc(c, box);
    return dot(c, c);
}

double dot(Vec3 a, Vec3 b)
{
    return a[X] * b[X] + a[Y] * b[Y] + a[Z] * b[Z];
}

double magnitude(Vec3 x)
{
    return __sqrt(dot(x, x));
}

double volume(CubicBox box)
{
    return box[X] * box[Y] * box[Z];
}

/*
#ifdef USE_ONEAPI
sycl::event distance_kernel(sycl::queue q, Vec3 a, Vec3 b, CubicBox box, const double& d)
{
    
}
sycl::event distance2_kernel(sycl::queue q,Vec3 a, Vec3 b, CubicBox box, const double& d)
{
    throw std::runtime_error("Not im");
}
sycl::event dot_kernel(sycl::queue q, Vec3 a, Vec3 b, const double& d)
{
    throw std::runtime_error("Not im");
}
sycl::event magnitude_kernel(sycl::queue q, Vec3 x, const double& m)
{
    throw std::runtime_error("Not im");
}

sycl::event pbc_kernel(sycl::queue q, Vec3 a, CubicBox box, const Vec3& v)
{
    
    //oneapi::mkl::vm::sub(q,)
    //throw std::runtime_error("Not im");
}

sycl::event volume_kernel(sycl::queue q, CubicBox box, const double& v)
{
    throw std::runtime_error("Not im");
}

#endif
*/