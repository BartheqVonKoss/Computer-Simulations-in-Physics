//
//  main.cpp
//  L1
//
//  Created by Bartłomiej Kos on 09/04/2018.
//  Copyright © 2018 Bartłomiej Kos. All rights reserved.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

void Ex1();
double radius(double x, double y);
void Ex2();
double radius(double x, double y)
{
    return sqrt(x * x + y * y);
}

void Ex1()                                // runge - kutta 4th order method
{
    long double h = 0.02;
    const int n = 100000;
    long double kx1, kx2, kx3, kx4, ky1, ky2, ky3, ky4;
    long double kvx1, kvx2, kvx3, kvx4, kvy1, kvy2, kvy3, kvy4;
    std::ofstream myData;
    myData.open("/Users/bartlomiejkos/Downloads/Programming/SymulacjeKomputeroweWFizyce/L1/data.csv", std::ios::out);
    std::vector<long double> x;
    std::vector<long double> y;
    std::vector<long double> vx;
    std::vector<long double> vy;
    std::vector<long double> EK;
    std::vector<long double> EP;
    x.push_back(0.5);
    y.push_back(0);
    vx.push_back(0);
    vy.push_back(1.63);
    for(int i=0; i<n; i++)                      // forces calculation
    {
        kx1 = vx[i];
        ky1 = vy[i];
        kvx1 = -x[i] / pow(radius(x[i], y[i]), 3);
        kvy1 = -y[i] / pow(radius(x[i], y[i]), 3);

        kx2 = vx[i] + kvx1 * h / 2;
        ky2 = vy[i] + kvy1 * h / 2;
        kvx2 = -(x[i] + h * kx1 / 2) / pow(radius(x[i] + kx1 * h / 2, y[i] + ky1 * h / 2), 3);
        kvy2 = -(y[i] + h * ky1 / 2) / pow(radius(x[i] + kx1 * h / 2, y[i] + ky1 * h / 2), 3);

        kx3 = vx[i] + kvx2 * h / 2;
        ky3 = vy[i] + kvy2 * h / 2;
        kvx3 = -(x[i] + h * kx2 / 2) / pow(radius(x[i] + kx2 * h / 2, y[i] + ky2 * h / 2), 3);
        kvy3 = -(y[i] + h * ky2 / 2) / pow(radius(x[i] + kx2 * h / 2, y[i] + ky2 * h / 2), 3);

        kx4 = vx[i] + kvx3 * h;
        ky4 = vy[i] + kvy3 * h;
        kvx4 = -(x[i] + h * kx3) / pow(radius(x[i] + kx3 * h, y[i] + ky3 * h), 3);
        kvy4 = -(y[i] + h * ky3) / pow(radius(x[i] + kx3 * h, y[i] + ky3 * h), 3);


        x.push_back(x[i] + h / 6 * (kx1 + 2 * kx2 + 2 * kx3 + kx4));                // data appending to determine positions and velocities calculated with forces
        y.push_back(y[i] + h / 6 * (ky1 + 2 * ky2 + 2 * ky3 + ky4));
        vx.push_back(vx[i] + h / 6 * (kvx1 + 2 * kvx2 + 2 * kvx3 + kvx4));
        vy.push_back(vy[i] + h / 6 * (kvy1 + 2 * kvy2 + 2 * kvy3 + kvy4));

        EP.push_back(-pow(sqrt(x[i] * x[i] + y[i] * y[i]), -1));
        EK.push_back(pow(sqrt(vx[i] * vx[i] + vy[i] * vy[i]), 2));
        std::cout << x[i] << '\t' << y[i] << '\t' << vx[i] << '\t' << vy[i] << '\t' << EP[i] << '\t' << EK[i] << '\t' << EK[i] + EP[i] << '\t' << i * h << std::endl; // printing & to file exporting
        myData << x[i] << ',' << y[i] << ',' << vx[i] << ',' << vy[i] << ',' << EP[i] << ',' << EK[i] << ',' << EK[i] + EP[i] << ',' << i * h << '\n';
    }

}

void Ex2()                                // nose oscilations on the plane with runge - kutta
{
    std::ofstream myData;
    myData.open("/Users/bartlomiejkos/Downloads/Programming/SymulacjeKomputeroweWFizyce/L1/2data.csv", std::ios::out);
    int n = 50000;
    double Q = 0.1;                       // to change for diffrent variant of nose oscilations
    std::vector<double> q;
    std::vector<double> p;
    std::vector<double> z;
    q.push_back(0);
    p.push_back(1);
    z.push_back(0);
    double h = 0.02;
    double kq1, kq2, kq3, kq4;
    double kp1, kp2, kp3, kp4;
    double kz1, kz2, kz3, kz4;

    for (int i = 0; i < n; i++)
    {
        kq1 = h * p[i];
        kp1 = h * (-q[i] - z[i] * p[i]);
        kz1 = h * ((p[i] * p[i] - 1) / Q);

        kq2 = h * p[i] + h * kq1 / 2;
        kp2 = h * (-q[i] - z[i] * p[i]) + h * kp1 / 2;
        kz2 = h * ((p[i] * p[i] - 1) / Q) + h * kz1 / 2;

        kq3 = h * p[i] + h * kq2 / 2;
        kp3 = h * (-q[i] - z[i] * p[i]) + h * kp2 / 2;
        kz3 = h * ((p[i] * p[i] - 1) / Q) + h * kz2 / 2;

        kq4 = h * p[i] + h * kq3;
        kp4 = h * (-q[i] - z[i] * p[i]) + h * kp3;
        kz4 = h * ((p[i] * p[i] - 1) / Q) + h * kz3;

        q.push_back(q[i] + h / 6 * (kq1 + 2 * kq2 + 2 * kq3 + kq4));
        p.push_back(p[i] + h / 6 * (kp1 + 2 * kp2 + 2 * kp3 + kp4));
        z.push_back(z[i] + h / 6 * (kz1 + 2 * kz2 + 2 * kz3 + kz4));
    }
    for (int i = 0; i < n; i++)
    {
        std::cout << p[i] << '\t' << q[i] << std::endl;
        myData << p[i] << ',' << q[i] << '\n';
    }
}

int main(int argc, const char * argv[])
{
    Ex1();
    Ex2();
    return 0;
}
