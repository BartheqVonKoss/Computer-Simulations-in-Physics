//
//  bivector.cpp
//  RDF_3rd
//
//  Created by Bartłomiej Kos on 29/05/2018.
//  Copyright © 2018 Bartłomiej Kos. All rights reserved.
//

#include "bivector.hpp"

BiVector::BiVector() : x_pos(0), y_pos(0) {};

BiVector::BiVector(double x, double y)
{
    this->x_pos = x;
    this->y_pos = y;
}

std::ostream &operator<<(std::ostream &F, const BiVector &in)
{
    F << "(" << in.x_pos << "," << in.y_pos << ")" << '\n';
    return F;
}

PParticle::PParticle(double x, double y) : position(x, y), velocity(), force() {}

std::ostream &operator<<(std::ostream &F, const PParticle &in)
{
    F << "(" << in.position.x_pos << "," << in.position.y_pos << ")" << "  (" << in.velocity.x_pos << "," << in.velocity.y_pos << ")  (" << in.force.x_pos << "," << in.force.y_pos << ")" << '\n';
    return F;
}
