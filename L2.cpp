//
//  main.cpp
//  L2
//
//  Created by Bartłomiej Kos on 07/04/2018.
//  Copyright © 2018 Bartłomiej Kos. All rights reserved.
//

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

double radius(double x, double y)
{
    return (sqrt(x * x + y * y));
}

double Fx(double x, double r)
{
    return (- x / (r * r * r));
}

double Fy(double y, double r)
{
    return (- y / (r * r * r));
}

void VerletAlgorithm()
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> vx;
    std::vector<double> vy;
    x.push_back(0.5);
    y.push_back(0);
    vx.push_back(0);
    vy.push_back(1.63);
    int n = 1000;
    double h = 0.02;
    double r = radius(x[0], y[0]);
    double ax = Fx(x[0], r);
    double ay = Fx(y[0], r);
    x.push_back(x[0] - h * vx[0] + h * h * ax / 2);
    y.push_back(y[0] - h * vy[0] + h * h * ay / 2);
    for(int i = 1; i < n; i++)
    {
        r = radius(x[i], y[i]);
        ax = Fx(x[i], r);
        ay = Fy(y[i], r);
        x.push_back(2 * x[i] - x[i-1] + h * h * ax);
        y.push_back(2 * y[i] - y[i-1] + h * h * ay);
        vx.push_back((x[i+1] - x[i-1]) / (2 * h));
        vy.push_back((y[i+1] - y[i-1]) / (2 * h));
    }
    for(int i = 0; i < n; i++)
    {
        std::cout << x[i] << '\t' << y[i] << '\t' << vx[i] << '\t' << vy[i] << '\t' << std::endl;
    }
}

void VerletVelocity()
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> vx;
    std::vector<double> vy;
    x.push_back(0.5);
    y.push_back(0);
    vx.push_back(0);
    vy.push_back(1.63);
    int n = 10000;
    double h = 0.002;
    for(int i = 0; i < n; i++)
    {
        x.push_back(x[i] + h * vx[i] + h * h * Fx(x[i], radius(x[i], y[i])));
        y.push_back(y[i] + h * vy[i] + h * h * Fy(y[i], radius(x[i], y[i])));
        vx.push_back(vx[i] + h * (Fx(x[i+1], radius(x[i+1], y[i+1])) + Fx(x[i], radius(x[i], y[i]))));
        vy.push_back(vy[i] + h * (Fy(y[i+1], radius(x[i+1], y[i+1])) + Fy(y[i], radius(x[i], y[i]))));
    }
    for(int i = 0; i < n; i++)
    {
        std::cout << x[i] << '\t' << y[i] << '\t' << vx[i] << '\t' << vy[i] << '\t' << std::endl;
    }
}

int main(int argc, const char * argv[])
{
    //VerletAlgorithm();
    VerletVelocity();
    return 0;
}
