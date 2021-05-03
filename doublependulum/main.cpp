#include <iostream>
#include <iterator>
#include <cmath>
#include <vector>
#include "miniwindow.h"
#include "vector2.h"


double const g = 9.81;
std::vector<double> k1(4);   //a Runge-Kuttához szükséges k változók
std::vector<double> k2(4);
std::vector<double> k3(4);
std::vector<double> k4(4);

class pendulum{
public:
    double m1 = 1.5753;    //kezdeti feltételek
    double m2 = 0.5753;
    double l1 = 0.5*200;
    double l2 = 0.5*200;
    double phi1 = M_PI/2-0.3;
    double phi2 = -M_PI/4;
    double dphi1 = 0;
    double dphi2 = 0;
    double h = 0.05;

    double ddphi1(double x1, double x2, double dx1, double dx2){       //a második deriváltakat kiértékelő függvények
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

    void RK4(){     //a jegyzőkönyvben leírt módon a Runge-Kutta
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

        if (phi1 > M_PI){phi1-=2*M_PI;}    //a fázistér miatt, hogyha átfordul az inga, akkor nem növeli tovább a phi értékét
        if (phi1 < -M_PI){phi1+=2*M_PI;}
        if (phi2 > M_PI){phi2-=2*M_PI;}
        if (phi2 < -M_PI){phi2+=2*M_PI;}
    }
};


int main() {
    int x = 0, y = 0;
    double mphi1 = 0, mphi2 = 0, mdphi1 = 0, mdphi2=0;  //az egér és a két inga által bezárt szögek, és szögsebességek
    bool a=false, b=false;                              //ezekben lesz tárolva, hogy le van e nyomva az egér

    pendulum pend;
    std::vector<std::vector<vector2d<double>>> data{{{pend.phi1,pend.dphi1},{pend.phi2,pend.dphi2}}};  //minden lefutáskor ebbe írja a szögeket ésa szögsebességeket
                                                                                                        //a fázistér miatt van rá szükség
    MainWindow wnd;
    wnd.window.eventDriven = false;
    wnd.mouseHandler([&](Mouse const& m)
                     {
                         x = m.x; y = m.y;
                         //a következő sorok az egér és az ingák közötti szöget számítják ki, és a szögsebességet határozzák meg,
                         //amivel elinduljon az inga ha elengedik
                         //több esetre azért van szükség, mivel az előjelek miatt anélkül nem tudná megkülönböztetni például a bal felső és a jobb alsó térfelet
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

                         //érzékeli az egér lenyomását, és az a-t vagy  b-t true-val teszi egyeblőé, attól függően, hogy melyik ingához van közlebb az egér
                         if(m.event == Mouse::LeftDown  ){
                             if (pow(320+pend.l1*sin(pend.phi1)-x,2)+pow(200+pend.l1*cos(pend.phi1)-y,2)<
                                pow(320+pend.l1*sin(pend.phi1)+pend.l2*sin(pend.phi2)-x,2)+pow(200+pend.l1*cos(pend.phi1)+pend.l2*cos(pend.phi2)-y,2)){
                                 a = true;
                             }
                             else{b = true;}}
                         //érzékeli az egér elengedését, visszaállítja az a-t és a b-t és beállítja a sebességet
                         else if(m.event == Mouse::LeftUp    ){
                             if(b){pend.dphi2 = -10*mdphi2;}
                             if(a){pend.dphi1 = -10*mdphi1;}
                             a=false;
                             b=false;
                         }
                     });

    wnd.resizeHandler([&](int w, int h, StateChange sc){ } );
    wnd.idleHandler([&]{

        //ha az a vagy a b true akkor azokat irányítja az egér, ha nem, futtatja az RK-t
        if (b){
            pend.phi2 = mphi2;
            pend.dphi1 = 0;
        }
        else if (a){
            pend.phi1 = mphi1;
            pend.dphi2 = 0;
        }
        else {pend.RK4();}
        data.push_back({{pend.phi1,pend.dphi1},{pend.phi2,pend.dphi2}});

    });
    wnd.exitHandler([&]{ });
    wnd.renderHandler( [&](SoftwareRenderer& r){

        r.forall_pixels([](auto, auto, auto){ return color(255, 255, 255); });

        //a két inga ábrázolása
        r.ellipse(int(320),int(200),5,5,[](auto){ return color(100, 100, 100); });
        r.line(320,200,int(320+pend.l1*sin(pend.phi1)),int(200+pend.l1*cos(pend.phi1)),[](auto){ return color(100, 100, 100); });
        r.ellipse(int(320+pend.l1*sin(pend.phi1)),int(200+pend.l1*cos(pend.phi1)),5,5,[](auto){ return color(0, 200, 100); });
        r.line(int(320+pend.l1*sin(pend.phi1)),int(200+pend.l1*cos(pend.phi1)),int(320+pend.l1*sin(pend.phi1)+pend.l2*sin(pend.phi2)),int(200+pend.l1*cos(pend.phi1)+pend.l2*cos(pend.phi2)),[](auto){ return color(100, 100, 100); });
        r.ellipse(int(320+pend.l1*sin(pend.phi1)+pend.l2*sin(pend.phi2)),int(200+pend.l1*cos(pend.phi1)+pend.l2*cos(pend.phi2)),5,5,[](auto){ return color(255, 0, 0); });


        //fázistér
        for (auto e : data){
            r.setpixel(900+(int)100*e[0].x,200+(int)100*e[0].y,color(0, 200, 100));
            r.setpixel(900+(int)100*e[1].x,200+(int)100*e[1].y,color(255, 0, 0));
        }

        //koordináta rendszer
        r.line(580,200,1220,200,[](auto){ return color(100, 100, 100); });
        r.line(900,10,900,390,[](auto){ return color(100, 100, 100); });
        r.line(585,205,585,195,[](auto){ return color(100, 100, 100); });
        r.line(1214,205,1214,195,[](auto){ return color(100, 100, 100); });
        r.line(580,0,580,440,[](auto){ return color(0, 0, 0); });

    });

    bool res = wnd.open(L"Double pendulum", {64, 64}, {1280, 480},
                        true, [&]{ return true; });

    return res ? 0 : -1;

}
