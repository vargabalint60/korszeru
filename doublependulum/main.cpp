#include <iostream>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <vector>
#include "miniwindow.h"


double const g = 9.81;
std::vector<double> k1(4);
std::vector<double> k2(4);
std::vector<double> k3(4);
std::vector<double> k4(4);


class pendulum{
public:
    double m1 = 1;
    double m2 = 2;
    double l1 = 100;
    double l2 = 100;
    double phi1 = 1;
    double phi2 = 1;
    double dphi1 = -0.3;
    double dphi2 = 0;
    double h = 0.05;



    double ddphi1(double x1, double x2, double dx1, double dx2){
        return (m2*l1*dx1*dx1*sin(x2-x1)*cos(x2-x1) + m2*g*sin(x2)*cos(x2-x1) + m2*l2*dx2*dx2*sin(x2-x1) -
                (m1+m2)*g*sin(x1))/((m1+m2)*l1 - l1*cos(x2-x1)*cos(x2-x1));}

    double ddphi2(double x1, double x2, double dx1, double dx2){
        return (-m2*l2*dx2*dx2*sin(x2-x1)*cos(x2-x1) + (m1+m2)*(g*sin(x1)*cos(x2-x1) - l1*dx1*dx1*sin(x2-x1) -
                g*sin(x2)))/((m1+m2)*l2-m2*l2*cos(x2-x1)*cos(x2-x1));}

    void euler() {
        dphi1 += h * ddphi1(phi1, phi2, dphi1, dphi2);
        dphi2 += h * ddphi2(phi1, phi2, dphi1, dphi2);
        phi1 += h * dphi1;
        phi2 += h * dphi2;
    }

    void RK(){
        k1[0] = dphi1;
        k1[1] = ddphi1(phi1,phi2,dphi1,dphi2);
        k1[2] = dphi2;
        k1[3] = ddphi2(phi1,phi2,dphi1,dphi2);

        k2[0] = dphi1 + h/2*k1[0];
        k2[1] = ddphi1(phi1+h/2*k1[0],phi2+h/2*k1[2],dphi1+h/2*k1[1],dphi2+h/2*k1[3]);
        k2[2] = dphi2 + h/2*k1[2];
        k2[3] = ddphi2(phi1+h/2*k1[0],phi2+h/2*k1[2],dphi1+h/2*k1[1],dphi2+h/2*k1[3]);

        k3[0] = dphi1 + h/2*k2[0];
        k3[1] = ddphi1(phi1+h/2*k2[0],phi2+h/2*k2[2],dphi1+h/2*k2[1],dphi2+h/2*k2[3]);
        k3[2] = dphi2 + h/2*k2[2];
        k3[3] = ddphi2(phi1+h/2*k2[0],phi2+h/2*k2[2],dphi1+h/2*k2[1],dphi2+h/2*k2[3]);

        k4[0] = dphi1 + h*k3[0];
        k4[1] = ddphi1(phi1+h*k3[0],phi2+h*k3[2],dphi1+h*k3[1],dphi2+h*k3[3]);
        k4[2] = dphi2 + h*k3[2];
        k4[3] = ddphi2(phi1+h*k3[0],phi2+h*k3[2],dphi1+h*k3[1],dphi2+h*k3[3]);

        phi1 += h/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        dphi1 += h/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
        phi2 += h/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
        dphi2 += h/6*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3]);
    }

};


int main() {
    pendulum pend;
    MainWindow wnd;
    wnd.window.eventDriven = false;
    wnd.mouseHandler([&](Mouse const& m){ });
    wnd.resizeHandler([&](int w, int h, StateChange sc){ } );
    wnd.idleHandler([&]{
        // pend.euler();
        pend.RK();
    });
    wnd.exitHandler([&]{ });
    wnd.renderHandler( [&](SoftwareRenderer& r){
        r.forall_pixels([](auto, auto, auto){ return color(255, 255, 255); });

        r.ellipse(int(320),int(200),5,5,[](auto){ return color(100, 100, 100); });

        r.ellipse(int(320+pend.l1*sin(pend.phi1)),int(200+pend.l1*cos(pend.phi1)),5,5,[](auto){ return color(100, 100, 100); });

        r.ellipse(int(320+pend.l1*sin(pend.phi1)+pend.l2*sin(pend.phi2)),int(200+pend.l1*cos(pend.phi1)+pend.l2*cos(pend.phi2)),5,5,[](auto){ return color(100, 100, 100); });
    });

    bool res = wnd.open(L"Window name", {64, 64}, {640, 480},
                        true, [&]{ return true; });
    return res ? 0 : -1;


}
