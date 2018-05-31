//
//  bivector.hpp
//  RDF_3rd
//
//  Created by Bartłomiej Kos on 29/05/2018.
//  Copyright © 2018 Bartłomiej Kos. All rights reserved.
//

#ifndef bivector_hpp
#define bivector_hpp

#include <stdio.h>


#include <stdio.h>
#include <iostream>
#include <cmath>

struct BiVector
{
public:
    double x_pos;
    double y_pos;
public:
    BiVector();
    BiVector(double x, double y);
    //    void operator+=(const BiVector &in);
    //    void operator-=(const BiVector &in);
    friend std::ostream &operator<<(std::ostream &F, const BiVector &in);

    BiVector operator+=(BiVector &in)
    {
        this->x_pos += in.x_pos;
        this->y_pos += in.y_pos;
        return *this;
    }
    BiVector operator-=(BiVector &in)
    {
        this->x_pos -= in.x_pos;
        this->y_pos -= in.y_pos;
        return *this;
    }


    BiVector &operator*=(double s)
    {
        this->x_pos *= s;
        this->y_pos *= s;
        return *this;
    }
    BiVector &operator+=(double s)
    {
        this->x_pos += s;
        this->y_pos += s;
        return *this;
    }
    BiVector &operator-=(double s)
    {
        this->x_pos -= s;
        this->y_pos -= s;
        return *this;
    }
    BiVector &operator/=(double s)
    {
        this->x_pos /= s;
        this->y_pos /= s;
        return *this;
    }


};

struct PParticle
{
public:
    BiVector position;
    BiVector velocity;
    BiVector force;
public:
    PParticle(double x, double y);
    friend std::ostream &operator<<(std::ostream &F, const PParticle &in);

};



#endif /* bivector_hpp */
