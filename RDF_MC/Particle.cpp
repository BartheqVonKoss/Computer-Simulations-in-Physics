//
//  Particle.cpp
//  RDF_MC
//
//  Created by Bartłomiej Kos on 07/06/2018.
//  Copyright © 2018 Bartłomiej Kos. All rights reserved.
//

#include "Particle.hpp"
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <fstream>
#include <vector>
#include <random>
#include <cstdlib>


std::ostream &operator<<(std::ostream &F, const Particle &in)
{
    F << in.x << ',' << in.y;
    return F;
}

double distance(Particle &in1, Particle &in2, double L)
{
    double xr_x = in2.x - in1.x;
    double xr_y = in2.y - in1.y;
    if(xr_x > 0.5 * L) xr_x -= L;
    else if(xr_x < -0.5 * L) xr_x += L;
    if(xr_y > 0.5 * L) xr_y -= L;
    else if(xr_y < -0.5 * L) xr_y += L;
    return sqrt(xr_x * xr_x + xr_y * xr_y);
}

void initialization(std::vector<Particle> &atoms, int num_atoms, double L)
{
    int AUX = floor(sqrt(num_atoms)) + 1;
    double dL = L / (AUX + 1);
    int j = 0;
    for(int i = 0; i < num_atoms; i++)
    {
        if((i % AUX == 0) && (i != 0)) j++;
        atoms.push_back(Particle((1 + i - j * AUX) * dL,(1 + j) * dL));
    }
}

void boundary_conditions(std::vector<Particle> &atoms, int atom_no, double L)
{
    if(atoms[atom_no].x > L) atoms[atom_no].x -= L * round(atoms[atom_no].x / L);
    else if(atoms[atom_no].x < 0) atoms[atom_no].x += L * round(atoms[atom_no].x / L);
    if(atoms[atom_no].y > L) atoms[atom_no].y -= L * round(atoms[atom_no].y / L);
    else if(atoms[atom_no].y < 0) atoms[atom_no].y += L * round(atoms[atom_no].y / L);
}

void calculate_distances(std::vector<Particle> &atoms, int num_atoms, std::vector<double> &distances, int atom_no, double L)
{
    for(int i = 0; i < num_atoms; i++)
    {
        if(i == atom_no) continue;
        else
        {
            distances[i] = distance(atoms[atom_no], atoms[i], L);
        }
    }
}

double calculate_energy(std::vector<Particle> &atoms, int num_atoms, int atom_no, double rc, std::vector<double> &distances)
{
    double energy = 0;
    for(int i = 0; i < num_atoms - 1; i++)
    {
        if(i == atom_no) continue;
        else
        {
            if(distances[i] <= rc)
            {
                energy += 4 * (std::pow(1 / distances[i], 12) - std::pow(1 / distances[i], 6));
            }
        }
    }
    return energy;
}

void set_test_position(Particle &in)
{
    double r = ((double) rand() / (RAND_MAX)) + 1;
    double lambda = 0.005;
    in.x = in.x + lambda * (r - 0.5);
    in.y = in.y + lambda * (r - 0.5);
}


void rdf(std::vector<Particle> atoms, int num_atoms, std::vector<Particle> &RDF_aggregated, double loops, int &sample_no)
{
    double iter = 0.01;
    double drw = 0.1;
    double dx, dy, r;
    double n = 0;
    double result, count_pairs, avg_pairs;
    while(n < loops - 1)
    {
        result = 0;
        count_pairs = 0;
        avg_pairs = 0;
        for(int i = 0; i < num_atoms - 1; i++)
        {
            for(int j = i + 1; j < num_atoms; j++)
            {
                dx = atoms[i].x - atoms[j].x;
                dy = atoms[i].y - atoms[j].y;

                r = sqrt(dx * dx + dy * dy);
                if(r > iter && r < iter + drw) count_pairs++;
            }
        }
        avg_pairs = (double)count_pairs / (double)num_atoms;
        result = avg_pairs / (2 * M_PI * iter * drw);
        RDF_aggregated[n].x = iter;
        RDF_aggregated[n].y += result;
        n++;
        iter += 0.01;
    }
    sample_no++;
}












