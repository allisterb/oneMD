

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

#include "Rdf.hh"

Rdf::Rdf() { }

Rdf::Rdf(int nbins, CubicBox &box, string outfile)
{
    this->nbins = nbins;
    this->g.resize(this->nbins, 0.0);
    this->binwidth = box[0] / (2.0 * this->nbins);
    this->n = 0;
    this->outfile = outfile;
}

void Rdf::sample(vector <Vec3> &x, CubicBox &box)
{
    this->n++;
    #pragma omp parallel
    {
        vector <double> g_thread(nbins, 0.0);
        #pragma omp for schedule(guided, CHUNKSIZE)
        for (int i = 0; i < x.size()-1; i++)
        {
            for (unsigned int j = i+1; j < x.size(); j++)
            {
                double d = distance(x.at(i), x.at(j), box);
                if (d < box[0]/2.0)
                {
                    int ig = d/this->binwidth;
                    g_thread.at(ig) += 2.0;
                }
            }
        }

        #pragma omp critical
        {
            for (int i = 0; i < nbins; i++)
            {
                this->g.at(i) += g_thread.at(i);
            }
        }

    }
    return;
}

void Rdf::normalize(int natoms, CubicBox &box)
{
    double norm_factor = 4.0/3.0 * M_PI * natoms * (natoms-1.0) * this->n * pow(this->binwidth,3) / volume(box);

    #pragma omp parallel for schedule(guided, CHUNKSIZE)
    for (int i = 0; i < this->nbins; i++)
    {
        double r = (double) i;
        double binvol = pow(r+1.0,3) - pow(r,3);
        g.at(i) /= (binvol * norm_factor);
    }

    return;
}

void Rdf::output()
{
    ofstream oFS(this->outfile.c_str());
    oFS << setprecision(6) << fixed;
    for (int i = 0; i < this->nbins; i++)
    {
        oFS << setw(20) << i*this->binwidth;
        oFS << setw(20) << this->g.at(i);
        oFS << endl;
    }
    oFS.close();
    return;
}
