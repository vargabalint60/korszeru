#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cmath>

double y = 0;
double dy = y*y+1;
double h = 0.0015;
double k1,k2,k3,k4;

double euler(){
    y += h*dy;
    dy = y*y+1;
    return y;
}

double runge(){
    k1 = dy;
    k2 = 1 + (y+h*k1/2)*(y+h*k1/2);
    k3 = 1 + (y+h*k2/2)*(y+h*k2/2);
    k4 = 1 + + (y+h*k3)*(y+h*k3);
    y += h/6*(k1 + 2*k2 + 2*k3 + k4);
    dy = y*y+1;
    return y;
}

int main() {
    std::vector<double> x(1000);
    std::vector<double> x2(1000);

    x[0] = 0;
    x2[0] = 0;
    std::generate(x2.begin()+1,x2.end(),euler);

    y = 0;
    dy = y*y+1;
    std::generate(x.begin()+1,x.end(),runge);

    y = 0;
    std::ofstream output("data.txt");
    if( output.is_open() )
    {
        for (int i=0; i<1000; i++){
            output << x2[i] << ' ' << x[i] << " " << tan(y) << std::endl;
            y+=h;
        }
    }
    else{ std::cout << "Could not open output file\n"; }

    return 0;
}
