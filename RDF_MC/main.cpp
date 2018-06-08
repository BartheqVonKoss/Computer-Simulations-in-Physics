//
//  main.cpp
//  RDF_MC
//
//  Created by Bartłomiej Kos on 03/06/2018.
//  Copyright © 2018 Bartłomiej Kos. All rights reserved.
//

#include "Particle.hpp"



int main(int argc, const char * argv[])
{


    
    
    const double T = 2;                            // Some initial data
//    const double kB = 1.38064852e-23;
//    const double sigma = 0.34;
//    const double epsilon = 120*kB;
//    const double mass = 6.63e-23;
//    const double dt = 0.000064;
    double rc = 4.8;
    const double L = 10;
    const double num_atoms = 90;                     // No of atoms
    const double num_steps = 500000;                 // No of steps
    int sample_no = 0;
    double energy_container[2];
    std::ofstream myData, rdfData;
    myData.open("/Users/bartlomiejkos/Documents/Programming/C++/Wprowadzenie do C++/Symulacj komputerowe w fizyce/RDF_MC/data09T2.csv", std::ios::out);
    rdfData.open("/Users/bartlomiejkos/Documents/Programming/C++/Wprowadzenie do C++/Symulacj komputerowe w fizyce/RDF_MC/RDFdata08T2.csv", std::ios::out);
    std::vector<Particle> atoms;
    std::vector<double> distances;
    distances.resize(num_atoms);
    initialization(atoms, num_atoms, L);

    double loops = rc / 0.01;
    std::vector<Particle> RDF_aggregated(loops);
    rdf(atoms, num_atoms, RDF_aggregated, loops, sample_no);
    double cp =0;
    double ca =0;
    int atom_no;
    double dE;
    for(int i = 0; i < num_atoms; i++) myData << atoms[i] << std::endl;                                 // Print initial location of particles
    for(int j = 0; j < num_steps; j++)
    {
        if(j > 40000 && j % 100 == 0) rdf(atoms, num_atoms, RDF_aggregated, loops, sample_no);
        for(int k = 0; k < num_atoms; k++)
        {
            atom_no = (int)(((double) rand() / (RAND_MAX))  * num_atoms);                                            // Choose a particle
            calculate_distances(atoms, num_atoms, distances, atom_no, L);                               // Measure distances between the particles
            energy_container[0] = calculate_energy(atoms, num_atoms, atom_no, rc, distances);           // Hold value of energy for the particle before moving it
            Particle last_position = atoms[atom_no];                                                    // Remember position of particle before randomly generating its new place
            set_test_position(atoms[atom_no]);                                                          // Move the particle
            boundary_conditions(atoms, atom_no, L);                                                     // Periodic boundary conditions
            calculate_distances(atoms, num_atoms, distances, atom_no, L);                               // Measure distances between all particles
            energy_container[1] = calculate_energy(atoms, num_atoms, atom_no, rc, distances);           // Hold value of energy for the particle after moving it
            dE = energy_container[1] - energy_container[0];
            if(dE <= 0 || ((double) rand() / (RAND_MAX))  < exp(-dE / T))                               // Metropolis criterion
            {
                energy_container[0] = energy_container[1];
                cp++;
            }
            else
            {
                atoms[atom_no] = last_position;
            }
            ca++;
        }
    }
    for(int i = 0; i < num_atoms; i++)
    {
        myData << atoms[i] << std::endl;                                                                // Print final location of particles
    }
    for(int i = 0; i < loops - 1; i++)
    {
        RDF_aggregated[i].y /= sample_no;                                                               // Calculation of aggregated RDF(r) / no_of_rdf
        rdfData << RDF_aggregated[i] << '\n';
    }
    myData.close();
    rdfData.close();
    std::cout << cp / ca << std::endl;                                                                  // Calculation of success rate in Metropolis criterion
    return 0;
}
