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
 * @brief Header for Vec3 class
 * @see Vec3.h
 */

#include "Vec3.hh"

Vec3::Vec3(){ }

Vec3::Vec3(double x, double y, double z)
{
    this->r[X] = x;
    this->r[Y] = y;
    this->r[Z] = z;
}

double& Vec3::operator[](int i)
{
    return r[i];
}

const double& Vec3::operator[](int i) const
{
    return r[i];
}

void Vec3::set(double x, double y, double z)
{
    this->r[X] = x;
    this->r[Y] = y;
    this->r[Z] = z;
    return;
}

Vec3 Vec3::operator-(Vec3 rhs)
{
    return (Vec3 (r[X] - rhs[X], r[Y] - rhs[Y], r[Z] - rhs[Z]));
}

void Vec3::operator-=(Vec3 rhs)
{
    r[X] -= rhs[X];
    r[Y] -= rhs[Y];
    r[Z] -= rhs[Z];
    return;
}

Vec3 Vec3::operator+(Vec3 rhs)
{
    return (Vec3 (r[X] + rhs[X], r[Y] + rhs[Y], r[Z] + rhs[Z]));
}

void Vec3::operator+=(Vec3 rhs)
{
    r[X] += rhs[X];
    r[Y] += rhs[Y];
    r[Z] += rhs[Z];
    return;
}

Vec3 Vec3::operator/(double rhs)
{
    return (Vec3 (r[X] / rhs, r[Y] / rhs, r[Z] / rhs));
}

void Vec3::operator/=(double rhs)
{
    r[X] /= rhs;
    r[Y] /= rhs;
    r[Z] /= rhs;
    return;
}

Vec3 operator*(Vec3 lhs, double rhs)
{
    return (Vec3 (lhs[X] * rhs, lhs[Y] * rhs, lhs[Z] * rhs));
}

Vec3 operator*(double lhs, Vec3 rhs)
{
    return (Vec3 (rhs[X] * lhs, rhs[Y] * lhs, rhs[Z] * lhs));
}

void Vec3::operator*=(double rhs)
{
    r[X] *= rhs;
    r[Y] *= rhs;
    r[Z] *= rhs;
    return;
}

void Vec3::operator=(double rhs)
{
    r[X] = rhs;
    r[Y] = rhs;
    r[Z] = rhs;
    return;
}

