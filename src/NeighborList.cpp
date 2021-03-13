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

#include "NeighborList.hh"
#ifdef USE_ONEAPI
#include "dpc_common.hpp"
#include <CL/sycl.hpp>
using namespace oneapi;
#endif

NeighborList::NeighborList():
list(),
listptr(nullptr),
rlist(0.0),
rlist2(0.0)
{}

NeighborList::NeighborList(int natoms, double _rlist) :
list(natoms, vector<int>(natoms)),
listptr(nullptr),
rlist(_rlist),
rlist2(_rlist * _rlist)
{
    listptr = &list[0][0];
}

void NeighborList::UpdateHostCPU(vector <Vec3> &x, CubicBox &box)
{
#ifndef USE_ONEAPI
    #pragma omp parallel for schedule(guided, CHUNKSIZE)
    for (int i = 0; i < this->list.size(); i++)
    {
        this->list.at(i).resize(0);
    }
    #pragma omp parallel for schedule(guided, CHUNKSIZE)
    for (int i = 0; i < x.size()-1; i++)
    {
        for (unsigned int j = i+1; j < x.size(); j++)
        {
            if (distance2(x.at(i), x.at(j), box) < rlist2)
            {
                this->list.at(i).push_back(j);
            }
        }
    }
#else
    std::for_each(dpl::execution::par_unseq, this->list.begin(), this->list.end(), 
        [](vector <int> &v){v.resize(0);});
    // Atoms are not double counted in the neighbor list. That is, when atom j
    // is on atom i's list, the opposite is not true.

#endif    
return;
}


int NeighborList::GetSize(int i)
{
    return list.at(i).size();
}

int NeighborList::GetNeighbor(int i, int j)
{
    return list.at(i).at(j);
}
