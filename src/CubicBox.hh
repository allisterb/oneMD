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

#ifndef CubicBox_H
#define CubicBox_H

#include "Vec3.hh"
#include "xdrfile.h"
#include <vector>
using namespace std;

/** @brief Box dimensions.
 * @details This is just a two dimensional array initialized to three
 * items in each dimension. To access the elements of the array use operator().
 * For example, to if the box is cubic and your have a CubicBox object named
 * mybox, to get the X dimension do mybox(0). If you it is truly a CubicBox
 * (not cubic) you can access elements with mybox(i,j).*/
class CubicBox {
private:
    array <float,3> box;
public:
    CubicBox();
    CubicBox(float x, float y, float z);
    float& operator[] (int i);
    const float& operator[] (int i) const;

};

double distance(Vec3 a, Vec3 b, CubicBox box);
double distance2(Vec3 a, Vec3 b, CubicBox box);
double dot(Vec3 a, Vec3 b);
double magnitude(Vec3 x);
Vec3 pbc(Vec3 a, CubicBox box);
double volume(CubicBox box);

#endif
