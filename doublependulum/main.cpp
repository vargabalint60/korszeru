#include <iostream>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <vector>
#include "miniwindow.h"
#include "vector2.h"


double const g = 9.81;
std::vector<double> k1(4);
std::vector<double> k2(4);
std::vector<double> k3(4);
std::vector<double> k4(4);
std::vector<double> k5(4);
std::vector<double> k6(4);


class pendulum{
public:
    double m1 = 1;
    double m2 = 1;
    double l1 = 90;
    double l2 = 90;
    double phi1 = -0.4;
    double phi2 = -0.4;
    double dphi1 = 0;
    double dphi2 = 0;
    double h = 0.05;
    double z1=phi1,z2=phi2,z3=dphi1,z4=dphi2;



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

    void RK4(){
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

        if (phi1 > M_PI){phi1-=2*M_PI;}
        if (phi1 < -M_PI){phi1+=2*M_PI;}
        if (phi2 > M_PI){phi2-=2*M_PI;}
        if (phi2 < -M_PI){phi2+=2*M_PI;}
    }

    void RK45(){
        k1[0] = h*dphi1;
        k1[1] = h*ddphi1(phi1,phi2,dphi1,dphi2);
        k1[2] = h*dphi2;
        k1[3] = h*ddphi2(phi1,phi2,dphi1,dphi2);

        k2[0] = h*(dphi1 + 1.0/4*k1[0]);
        k2[1] = h*ddphi1(phi1+1.0/4*k1[0],phi2+1.0/4*k1[2],dphi1+1.0/4*k1[1],dphi2+1.0/4*k1[3]);
        k2[2] = h*(dphi2 + 1.0/4*k1[2]);
        k2[3] = h*ddphi2(phi1+1.0/4*k1[0],phi2+1.0/4*k1[2],dphi1+1.0/4*k1[1],dphi2+1.0/4*k1[3]);

        k3[0] = h*(dphi1 + 13.0/32.0*k1[0] + 9.0/32.0*k2[0]);
        k3[1] = h*ddphi1(phi1 + 13.0/32.0*k1[0] + 9.0/32.0*k2[0],phi2+13.0/32.0*k1[2] + 9.0/32.0*k2[2],
                         dphi1+13.0/32.0*k1[1] + 9.0/32.0*k2[1],dphi2+13.0/32.0*k1[3] + 9.0/32.0*k2[3]);
        k3[2] = h*(dphi2 + 13.0/32.0*k1[2] + 9.0/32.0*k2[2]);
        k3[3] = h*ddphi2(phi1 + 13.0/32.0*k1[0] + 9.0/32.0*k2[0],phi2+13.0/32.0*k1[2] + 9.0/32.0*k2[2],
                         dphi1+13.0/32.0*k1[1] + 9.0/32.0*k2[1],dphi2+13.0/32.0*k1[3] + 9.0/32.0*k2[3]);

        k4[0] = h*(dphi1 + 1932.0/2197.0*k1[0] - 7200.0/2197.0*k2[0] + 7296.0/2197.0*k3[0] );
        k4[1] = h*ddphi1(phi1+1932.0/2197.0*k1[0] - 7200.0/2197.0*k2[0] + 7296.0/2197.0*k3[0],phi2+1932.0/2197.0*k1[2] - 7200.0/2197.0*k2[2] + 7296.0/2197.0*k3[2],
                         dphi1+1932.0/2197.0*k1[1] - 7200.0/2197.0*k2[1] + 7296.0/2197.0*k3[1],dphi2+1932.0/2197.0*k1[3] - 7200.0/2197.0*k2[3] + 7296.0/2197.0*k3[3]);
        k4[2] = h*(dphi2 + 1932.0/2197.0*k1[2] - 7200.0/2197.0*k2[2] + 7296.0/2197.0*k3[2]);
        k4[3] = h*ddphi2(phi1+1932.0/2197.0*k1[0] - 7200.0/2197.0*k2[0] + 7296.0/2197.0*k3[0],phi2+1932.0/2197.0*k1[2] - 7200.0/2197.0*k2[2] + 7296.0/2197.0*k3[2],
                         dphi1+1932.0/2197.0*k1[1] - 7200.0/2197.0*k2[1] + 7296.0/2197.0*k3[1],dphi2+1932.0/2197.0*k1[3] - 7200.0/2197.0*k2[3] + 7296.0/2197.0*k3[3]);

        k5[0] = h*(dphi1 + 439.0/216.0*k1[0] - 8*k2[0] + 3680.0/513*k3[0] - 845.0/4104.0*k4[0]);
        k5[1] = h*ddphi1(phi1+ 439.0/216.0*k1[0] - 8*k2[0] + 3680.0/513*k3[0] - 845.0/4104.0*k4[0],phi2+ 439.0/216.0*k1[2] - 8*k2[2] + 3680.0/513*k3[2] - 845.0/4104.0*k4[2],
                         dphi1+ 439.0/216.0*k1[1] - 8*k2[1] + 3680.0/513*k3[1] - 845.0/4104.0*k4[1],dphi2+ 439.0/216.0*k1[3] - 8*k2[3] + 3680.0/513*k3[3] - 845.0/4104.0*k4[3]);
        k5[2] = h*(dphi2 + 439.0/216.0*k1[2] - 8*k2[2] + 3680.0/513*k3[2] - 845.0/4104.0*k4[2]);
        k5[3] = h*ddphi2(phi1+h*k3[0],phi2+h*k3[2],dphi1+h*k3[1],dphi2+h*k3[3]);

        k6[0] = h*(dphi1 - 8.0/27.0*k1[0] + 2*k2[0] - 3644.0/2565.0*k3[0] + 1859.0/4104.0*k4[0] - 11.0/40.0*k5[0]);
        k6[1] = h*ddphi1(phi1- 8.0/27.0*k1[0] + 2*k2[0] - 3644.0/2565.0*k3[0] + 1859.0/4104.0*k4[0] - 11.0/40.0*k5[0],phi2- 8.0/27.0*k1[2] + 2*k2[2] - 3644.0/2565.0*k3[2] + 1859.0/4104.0*k4[2] - 11.0/40.0*k5[2],
                       dphi1- 8.0/27.0*k1[1] + 2*k2[1] - 3644.0/2565.0*k3[1] + 1859.0/4104.0*k4[1] - 11.0/40.0*k5[1],dphi2- 8.0/27.0*k1[3] + 2*k2[3] - 3644.0/2565.0*k3[3] + 1859.0/4104.0*k4[3] - 11.0/40.0*k5[3]);
        k6[2] = h*(dphi2 - 8.0/27.0*k1[2] + 2*k2[2] - 3644.0/2565.0*k3[2] + 1859.0/4104.0*k4[2] - 11.0/40.0*k5[2]);
        k6[3] = h*ddphi2(phi1- 8.0/27.0*k1[0] + 2*k2[0] - 3644.0/2565.0*k3[0] + 1859.0/4104.0*k4[0] - 11.0/40.0*k5[0],phi2- 8.0/27.0*k1[2] + 2*k2[2] - 3644.0/2565.0*k3[2] + 1859.0/4104.0*k4[2] - 11.0/40.0*k5[2],
                       dphi1- 8.0/27.0*k1[1] + 2*k2[1] - 3644.0/2565.0*k3[1] + 1859.0/4104.0*k4[1] - 11.0/40.0*k5[1],dphi2- 8.0/27.0*k1[3] + 2*k2[3] - 3644.0/2565.0*k3[3] + 1859.0/4104.0*k4[3] - 11.0/40.0*k5[3]);

        z1 += 25.0/216.0*k1[0] + 1408.0/2565.0*k3[0] + 2197.0/4101.0*k4[0] - 1.0/5.0*k5[0];
        z2 += 25.0/216.0*k1[1] + 1408.0/2565.0*k3[1] + 2197.0/4101.0*k4[1] - 1.0/5.0*k5[1];
        z3 += 25.0/216.0*k1[2] + 1408.0/2565.0*k3[2] + 2197.0/4101.0*k4[2] - 1.0/5.0*k5[2];
        z4 += 25.0/216.0*k1[3] + 1408.0/2565.0*k3[3] + 2197.0/4101.0*k4[3] - 1.0/5.0*k5[3];

        phi1 += 16.0/135.0*k1[0] + 6656.0/12825.0*k3[0] + 28561.0/56430.0*k4[0] - 9.0/50.0*k5[0] + 2.0/55.0*k6[0];
        dphi1 += 16.0/135.0*k1[1] + 6656.0/12825.0*k3[1] + 28561.0/56430.0*k4[1] - 9.0/50.0*k5[1] + 2.0/55.0*k6[1];
        phi2 += 16.0/135.0*k1[2] + 6656.0/12825.0*k3[2] + 28561.0/56430.0*k4[2] - 9.0/50.0*k5[2] + 2.0/55.0*k6[2];
        dphi2 += 16.0/135.0*k1[3] + 6656.0/12825.0*k3[3] + 28561.0/56430.0*k4[3] - 9.0/50.0*k5[3] + 2.0/55.0*k6[3];

        h *= pow(1/(2*(fabs(phi1-z1) + fabs(phi2-z3) + fabs(dphi1-z2) + fabs(dphi2-z4))),1.0/4);

        //std::cout << pow(1/(2*(fabs(phi1-z1) + fabs(phi2-z3) + fabs(dphi1-z2) + fabs(dphi2-z4))),1.0/4) << std::endl;
    }

};


int main() {
    int x = 0, y = 0;
    double mphi1 = 0, mphi2 = 0, mdphi1 = 0, mphi3 = 0, mphi4 = 0, mdphi2=0;
    bool a=false, b=false;
    vector2d<double> state;

    pendulum pend;
    std::vector<std::vector<vector2d<double>>> data{{{pend.phi1,pend.dphi1},{pend.phi2,pend.dphi2}}};
    std::vector<vector2d<double>> data2;
    MainWindow wnd;
    wnd.window.eventDriven = false;
    wnd.mouseHandler([&](Mouse const& m)
                     {
                         x = m.x; y = m.y;
                         mphi3 = mphi1;
                         mphi4 = mphi2;
                         if (y < 200){
                             if (x > 320){mphi1 = atan((double)(x-320)/(y-200))+M_PI;}
                             if (x < 320){mphi1 = atan((double)(x-320)/(y-200))-M_PI;}
                         }
                         else{mphi1 = atan((double)(x-320)/(y-200));}
                         mdphi1 = pend.phi1-mphi1;

                         if (y < 200+pend.l1*cos(pend.phi1)){
                             if (x > 320+pend.l1*sin(pend.phi1)){mphi2 = atan((x-320-pend.l1*sin(pend.phi1))/(y-200-pend.l1*cos(pend.phi1)))+M_PI;}
                             if (x < 320+pend.l1*sin(pend.phi1)){mphi2 = atan((x-320-pend.l1*sin(pend.phi1))/(y-200-pend.l1*cos(pend.phi1)))-M_PI;}
                         }
                         else{mphi2 = atan((x-320-pend.l1*sin(pend.phi1))/(y-200-pend.l1*cos(pend.phi1)));}
                         mdphi2 = pend.phi2-mphi2;


                         if(m.event == Mouse::LeftDown  ){
                             if (pow(320+pend.l1*sin(pend.phi1)-x,2)+pow(200+pend.l1*cos(pend.phi1)-y,2)<
                                pow(320+pend.l1*sin(pend.phi1)+pend.l2*sin(pend.phi2)-x,2)+pow(200+pend.l1*cos(pend.phi1)+pend.l2*cos(pend.phi2)-y,2)){
                                 a = true;
                             }
                             else{b = true;}}
                         else if(m.event == Mouse::LeftUp    ){
                             if(b){pend.dphi2 = -10*mdphi2;}
                             if(a){pend.dphi1 = -10*mdphi1;}
                             a=false;
                             b=false;

                             //std::cout << mdphi1 << " " << mdphi2 << std::endl;
                         }


                     });

    wnd.resizeHandler([&](int w, int h, StateChange sc){ } );
    wnd.idleHandler([&]{
        // pend.euler();
        // pend.RK45();

        if (b){
            pend.phi2 = mphi2;
            pend.dphi1 = 0;
        }
        else if (a){
            pend.phi1 = mphi1;
            pend.dphi2 = 0;
        }
        else {pend.RK4();}
        //std::cout << pend.h << std::endl;

        //state = {pend.phi1,pend.dphi1};
        data.push_back({{pend.phi1,pend.dphi1},{pend.phi2,pend.dphi2}});

    });
    wnd.exitHandler([&]{ });
    wnd.renderHandler( [&](SoftwareRenderer& r){
        r.forall_pixels([](auto, auto, auto){ return color(255, 255, 255); });

        r.ellipse(int(320),int(200),5,5,[](auto){ return color(100, 100, 100); });
        r.line(320,200,int(320+pend.l1*sin(pend.phi1)),int(200+pend.l1*cos(pend.phi1)),[](auto){ return color(100, 100, 100); });
        r.ellipse(int(320+pend.l1*sin(pend.phi1)),int(200+pend.l1*cos(pend.phi1)),5,5,[](auto){ return color(0, 200, 100); });
        r.line(int(320+pend.l1*sin(pend.phi1)),int(200+pend.l1*cos(pend.phi1)),int(320+pend.l1*sin(pend.phi1)+pend.l2*sin(pend.phi2)),int(200+pend.l1*cos(pend.phi1)+pend.l2*cos(pend.phi2)),[](auto){ return color(100, 100, 100); });
        r.ellipse(int(320+pend.l1*sin(pend.phi1)+pend.l2*sin(pend.phi2)),int(200+pend.l1*cos(pend.phi1)+pend.l2*cos(pend.phi2)),5,5,[](auto){ return color(255, 0, 0); });

        for (auto e : data){
            r.setpixel(900+(int)100*e[0].x,200+(int)100*e[0].y,color(0, 200, 100));
            r.setpixel(900+(int)100*e[1].x,200+(int)100*e[1].y,color(255, 0, 0));
        }

        r.line(580,200,1220,200,[](auto){ return color(0, 0, 0); });
        r.line(900,10,900,390,[](auto){ return color(0, 0, 0); });
        r.line(585,205,585,195,[](auto){ return color(0, 0, 0); });
        r.line(1214,205,1214,195,[](auto){ return color(0, 0, 0); });
    });

    bool res = wnd.open(L"Double pendulum", {64, 64}, {1280, 480},
                        true, [&]{ return true; });

    return res ? 0 : -1;

}
