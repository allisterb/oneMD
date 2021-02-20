/*
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

/**
 * @file
 * @author James W. Barnett jbarnet4@tulane.edu
 * @date December 5, 2014
 * @brief Header for CubicBox class
 * @see CubicBox.h
 */

#include <fstream>
#include <math.h>
#include <string>
#include <vector>

#include "CubicBox.hh"

CubicBox::CubicBox() { }

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

Vector pbc(Vector a, CubicBox box)
{

    a[Z] -= box[Z] * nearbyint(a[Z] / box[Z]);
    a[Y] -= box[Y] * nearbyint(a[Y] / box[Y]);;
    a[X] -= box[X] * nearbyint(a[X] / box[X]);;
    return a;
}

double distance(Vector a, Vector b, CubicBox box)
{
    return sqrt(distance2(a, b, box));
}

double distance2(Vector a, Vector b, CubicBox box)
{
    Vector c = a - b;
    c = pbc(c, box);
    return dot(c, c);
}

double dot(Vector a, Vector b)
{
    return a[X] * b[X] + a[Y] * b[Y] + a[Z] * b[Z];
}

double magnitude(Vector x)
{
    return sqrt(dot(x, x));
}

double volume(CubicBox box)
{
    return box[X] * box[Y] * box[Z];
}
