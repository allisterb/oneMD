/* Based on libgmxcpp by James W. Barnett: https://github.com/wesbarnett/libgmxcpp
 * libgmxcpp
 * Copyright (C) 2015 James W. Barnett <jbarnet4@tulane.edu>
 s*
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

#include <fstream>
#include <string>
#include <vector> 

#include "common.hh"
#include "Vec3.hh"
#include "xdrfile.h"

using namespace std;

#ifndef USE_ONEAPI
class CubicBox {
    public:
        CubicBox();
        CubicBox(float x, float y, float z);
        float& operator[] (int i);
        const float& operator[] (int i) const;
        array <float,3> box;
};
#else
using CubicBox = sycl::double3;
#endif

SYCL_LINK Vec3 pbc(Vec3 a, CubicBox box);
SYCL_LINK double distance(Vec3 a, Vec3 b, CubicBox box);
SYCL_LINK double distance2(Vec3 a, Vec3 b, CubicBox box);
SYCL_LINK double dot(Vec3 a, Vec3 b);
SYCL_LINK double magnitude(Vec3 x);
SYCL_LINK double volume(CubicBox box);