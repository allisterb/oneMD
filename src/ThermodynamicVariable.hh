/* Based on https://github.com/wesbarnett/lennardjones" by James W. Barnett
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

#include <math.h>
#include <vector>

#include "Vec3.hh"

class ThermodynamicVariable {
    private:
        std::vector <double> all;
        double avg;
        double error;
    public:
        ThermodynamicVariable();
        double GetAvg();
        double GetError();
        void ErrorAnalysis(int nblocks);
        void Normalize();
        void Sample(double value);
};