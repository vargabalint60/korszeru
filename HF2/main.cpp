#include <iostream>
#include "vector2.h"

int main()
{
    vector2d<double> v1;
    std::cout << "Enter vector 1:" << std::endl;
    std::cin >> v1;
    std::cout << "v1 = " << v1 << std::endl;

    vector2d<double> v2;
    std::cout << "Enter vector 2:" << std::endl;
    std::cin >> v2;
    std::cout << "v2 = " << v2 << std::endl;
    vector2d<double> v3 = v2 - v1;
    std::cout << "v3 = v2-v1 = " << v3 << std::endl;

    vector2d<double> v4;
    std::cout << "v4 = " << v4 << std::endl;

    v3 += v2;
    std::cout << "v3 += v2 = " << v3 << std::endl;

    v3 -= v2;
    std::cout << "v3 -= v2 = " << v3 << std::endl;

    v3 /= 3;
    std::cout << "v3 /= 3 = " << v3 << std::endl;

    v3 *= 3;
    std::cout << "v3 *= 3 = " << v3 << std::endl;

    vector2d<double> v5 = 3*v3;
    std::cout << "v5 = 3*v3 = " << v5 << std::endl;

    v5 = v3*3;
    std::cout << "v5 = v3*3 = " << v5 << std::endl;

    double a = dot(v3,v2);
    std::cout << "dot(v3,v2) = " << a << std::endl;

    double l = length(v3);
    std::cout << "length(v3) = " << l << std::endl;

    double sql = sqlength(v3);
    std::cout << "sqlength(v3) = " << sql << std::endl;

    v4 = normalize(v3);
    std::cout << "v4 = normalize(v3) = " << v4 << std::endl;

    bool eq = v3 == v4;
    std::cout << "v3 == v4 = " << eq << std::endl;

    bool eq2 = v3 == v3;
    std::cout << "v3 == v3 = " << eq2 << std::endl;

    bool neq = v3 != v4;
    std::cout << "v3 != v4 = " << neq << std::endl;

    bool neq2 = v3 != v3;
    std::cout << "v3 != v3 = " << neq2 << std::endl;

    return 0;
}
