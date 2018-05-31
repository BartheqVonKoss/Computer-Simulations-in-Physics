//
//  main.cpp
//  RDF_3rd
//
//  Created by Bartłomiej Kos on 29/05/2018.
//  Copyright © 2018 Bartłomiej Kos. All rights reserved.
//

#include <iostream>
#include "bivector.hpp"
#include <vector>
#include <math.h>
#include <fstream>

void boundary_conditions(std::vector<PParticle> &atoms, int No, double L)
{

    if(atoms[No].position.x_pos > L) atoms[No].position.x_pos -= L * round(atoms[No].position.x_pos / L);
    else if(atoms[No].position.x_pos < 0) atoms[No].position.x_pos += L * round(atoms[No].position.x_pos / L);
    if(atoms[No].position.y_pos > L) atoms[No].position.y_pos -= L * round(atoms[No].position.y_pos / L);
    else if(atoms[No].position.y_pos < 0) atoms[No].position.y_pos += L * round(atoms[No].position.y_pos / L);

}

void initialization(std::vector<PParticle> &atoms, int num_atoms, double L)
{
    int AUX = floor(sqrt(num_atoms)) + 1;
    double dL = L / (AUX + 1);
    int j = 0;
    for(int i = 0; i < num_atoms - 1; i++)
    {
        if((i % AUX == 0) && (i != 0)) j++;
        atoms.push_back(PParticle((1 + i - j * AUX) * dL,(1 + j) * dL));
    }
}

void calculation_of_forces(std::vector<PParticle> &atoms, int num_atoms, double L)
{
    for(int i = 0; i < num_atoms; i++)
    {
        atoms[i].force.x_pos = 0;
        atoms[i].force.y_pos = 0;
    }
    double xr_x = 0, xr_y = 0, r2 = 0, ff = 0;
    for(int i = 0; i < num_atoms - 1; i++)
    {
        for(int j = i + 1; j < num_atoms; j++)
        {
            xr_x = atoms[i].position.x_pos - atoms[j].position.x_pos;
            xr_y = atoms[i].position.y_pos - atoms[j].position.y_pos;

            if(xr_x > 0.5 * L) xr_x -= L;
            else if(xr_x < -0.5 * L) xr_x += L;
            if(xr_y > 0.5 * L) xr_y -= L;
            else if(xr_y < -0.5 * L) xr_y += L;

            r2 = xr_x * xr_x + xr_y * xr_y;

            ff = std::pow(r2, -7) - 0.5 * std::pow(r2, -4);

            atoms[i].force.x_pos += ff * xr_x;
            atoms[i].force.y_pos += ff * xr_y;
            atoms[j].force.x_pos -= ff * xr_x;
            atoms[j].force.y_pos -= ff * xr_y;
        }
    }
}

void integration_of_motion(std::vector<PParticle> &atoms, int num_atoms, double dt, double L)
{
    for(int i = 0; i < num_atoms; i++)
    {

        atoms[i].position.x_pos += atoms[i].velocity.x_pos + 0.5 * atoms[i].force.x_pos * dt * dt;
        atoms[i].position.y_pos += atoms[i].velocity.y_pos + 0.5 * atoms[i].force.y_pos * dt * dt;

        boundary_conditions(atoms, i, L);

        atoms[i].velocity.x_pos += 0.5 * dt * atoms[i].force.x_pos;
        atoms[i].velocity.y_pos += 0.5 * dt * atoms[i].force.y_pos;

        calculation_of_forces(atoms, num_atoms, L);

        atoms[i].velocity.x_pos += 0.5 * dt * atoms[i].force.x_pos;
        atoms[i].velocity.y_pos += 0.5 * dt * atoms[i].force.y_pos;
    }
}

double RDF(std::vector<PParticle> &atoms, int num_atoms, double rw)
{
    const double sigma = 0.34;
    double V = 2 * M_PI * rw;
    int count_pairs = 0;
    double dx, dy, r2, r;
    double result;
    double drw = 0.1;
    double avg_pairs = 0;
    for(int i = 0; i < num_atoms - 1; i++)
    {
        for(int j = i + 1; j < num_atoms; j++)
        {
            dx = atoms[i].position.x_pos - atoms[j].position.x_pos;
            dy = atoms[i].position.y_pos - atoms[j].position.y_pos;

            r2 = dx * dx + dy * dy;
            r = sqrt(r2);
            if(r > rw && r < rw + drw) count_pairs++;
        }
    }
    avg_pairs = (double)count_pairs / (double)num_atoms;
    result = avg_pairs / (2 * M_PI * rw * drw);
    //result = avg_pairs / ((num_atoms * (num_atoms - 1)) * 2 * M_PI * rw * drw);
    return result;
}




int main(int argc, const char * argv[])
{
    const double kB = 1.38064852e-23;
    const double sigma = 0.34;
    const double epsilon = 120*kB;
    const double mass = 6.63e-23;
    const double dt = 0.000064;
    const double L = 10;
    const double num_atoms = 90;
    const double num_steps = 800;

    std::ofstream myData, rdfData;
    myData.open("/Users/bartlomiejkos/Documents/Programming/C++/Wprowadzenie do C++/Symulacj komputerowe w fizyce/RDF_3rd/data09.csv", std::ios::out);
    rdfData.open("/Users/bartlomiejkos/Documents/Programming/C++/Wprowadzenie do C++/Symulacj komputerowe w fizyce/RDF_3rd/RDFdata09.csv", std::ios::out);

    std::vector<PParticle> atoms;
    initialization(atoms, num_atoms, L);

    double rw = 0.1;
    while(rw < 4)
    {
        //std::cout << "liczba par (r) " << rw << " " << RDF(atoms, num_atoms, rw) << std::endl;
        rdfData << rw << ',' << RDF(atoms, num_atoms, rw) << '\n';
        rw += 0.01;
    }
    for(int i = 0; i < num_atoms; i++)
        myData << atoms[i].position.x_pos << ',' << atoms[i].position.y_pos<< '\n';

    calculation_of_forces(atoms, num_atoms, L);
    int j = 0;
    while(j < num_steps)
    {
        calculation_of_forces(atoms, num_atoms, L);
        integration_of_motion(atoms, num_atoms, dt, L);
        for(int i = 0; i < num_atoms; i++) std::cout << atoms[i].position.x_pos << ',' << atoms[i].position.y_pos<< '\n';
        j++;
    }
    std::cout << std::endl;
    rw = 0.1;
    double w = 0;
    while(rw < 4)
    {
        //std::cout << "liczba par (r) " << rw << " " << RDF(atoms, num_atoms, rw) << std::endl;
        rdfData << rw << ',' << RDF(atoms, num_atoms, rw) << '\n';

        rw += 0.01;
        w++;
        //std::cout << w << std::endl;
    }
    for(int i = 0; i < num_atoms; i++)
    {
        myData << atoms[i].position.x_pos << ',' << atoms[i].position.y_pos<< '\n';
        std::cout << atoms[i].position.x_pos << ',' << atoms[i].position.y_pos<< '\n';
    }
    return 0;
}
