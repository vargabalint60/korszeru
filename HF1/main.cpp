#include <iostream>
#include <cmath>

template <typename F, typename dF, typename T, typename R>
double Newton(F f, dF df, T x, R run){
    double x2 = x - f(x)/df(x) ;
    double x3;
    do {
        x3 = x2;
        x2 = x2 - f(x2)/df(x2);
    } while (run(x2,x3));
    return x2;
}

bool run(double x2, double x3){
    if (fabs(x3-x2) > 1e-5){return true;}
    else {return false;}
}

int main() {
    double root = Newton([](double x){ return x*x - 612.0; },
           [](double x){ return 2.0*x; }, 10.0, run);
    std::cout << root << std::endl;
    return 0;
}
