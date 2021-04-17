#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <iterator>

class fx{
    double F;
    double q;
    double a;
    double h;

public:
    fx(double F,double q,double a,double h) : F{F}, q{q}, a{a}, h{h}{}
    fx() : F{900}, q{1.8}, a{200}, h{35}{}

    double eval(double x){
        return F/q*(cosh(q/F*x)-cosh(a*q/(2*F)))+h;
    }

    double deriv(double x, double d){
        return (256*(eval(x+d) - eval(x-d)) - 40*(eval(x+2*d)-eval(x-2*d)) + eval(x+4*d)-eval(x-4*d))/(360*d);
    }

    double mid(int n){
        std::vector<double> l(n);
        std::vector<double> x(n);
        double dx = a/n;
        double x0=-a/2;
        auto dl = [&](auto xc){return sqrt(1+pow(deriv(xc,20),2))*dx;};

        std::generate(x.begin(),x.end(),[&]{return (x0+=dx)-0.5*dx;});
        std::transform(x.begin(),x.end(),l.begin(),dl);
        return std::accumulate(l.begin(),l.end(),(double)0);
    }

    double trap(int n) {
        std::vector<double> l(n + 1);
        std::vector<double> x(n + 1);
        double dx = a / n;
        double x0 = -a / 2;
        auto dl = [&](auto xc) { return sqrt(1 + pow(deriv(xc,20), 2)) * dx/2; };

        std::generate(x.begin(), x.end(), [&] { return (x0 += dx) - dx; });
        std::transform(x.begin(), x.end(), l.begin(), dl);
        std::transform(l.begin()+1,l.end()-1, l.begin(),[](auto lc){return 2*lc;});
        return std::accumulate(l.begin(),l.end(),(double)0);
    }

    double simpson(int n){
        return (2*mid(n)+trap(n))/3;
    }
};



int main() {
    fx f;
    int n = 5;
    std::vector<double> data = {f.mid(n),f.trap(n),f.simpson(n)};

    std::ofstream output("data.txt");
    if( output.is_open() )
    {
        output << "M_" << n <<" = " << data[0] << ", " << "T_" << n <<" = " << data[1] << ", " << "S_" << n <<" = " << data[2];
    }
    else{ std::cout << "Could not open output file\n"; }

    return 0;
}
