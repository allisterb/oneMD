/*
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


#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#ifdef USE_ONEAPI
// oneDPL headers should be included before standard headers
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
#include <oneapi/dpl/iterator>
#endif

#include "Vector.hh"
#include "CubicBox.hh"

class NeighborList {
    public:
        vector <vector<int>> list;
        double rlist;
        double rlist2;
        void* listptr;
        NeighborList();
        NeighborList(int natoms, double rlist);
        int GetNeighbor(int i, int j);
        int GetSize(int i);
        void Init(int natoms, double rlist);
        void UpdateHostCPU(vector <Vector> &x, CubicBox &box);
};

#endif
