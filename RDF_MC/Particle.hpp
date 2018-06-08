//
//  Particle.hpp
//  RDF_MC
//
//  Created by Bartłomiej Kos on 07/06/2018.
//  Copyright © 2018 Bartłomiej Kos. All rights reserved.
//

#ifndef Particle_hpp
#define Particle_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <fstream>
#include <vector>



class Particle
{
public:
    double x;
    double y;
public:
    Particle() : x(0), y(0) {};
    Particle(double xx, double yy) : x(xx), y(yy) {};
    Particle(const Particle &in) : x(in.x), y(in.y) {};
    friend std::ostream &operator<<(std::ostream &F, const Particle &in);
    double distance(Particle &in1, Particle &in2, double L);
};
void boundary_conditions(std::vector<Particle> &atoms, int atom_no, double L);
void calculate_distances(std::vector<Particle> &atoms, int num_atoms, std::vector<double> &distances, int atom_no, double L);
double calculate_energy(std::vector<Particle> &atoms, int num_atoms, int atom_no, double rc, std::vector<double> &distances);
void set_test_position(Particle &in);
void initialization(std::vector<Particle> &atoms, int num_atoms, double L);
void rdf(std::vector<Particle> atoms, int num_atoms, std::vector<Particle> &RDF_aggregated, double loops, int &sample_no);

#endif /* Particle_hpp */
