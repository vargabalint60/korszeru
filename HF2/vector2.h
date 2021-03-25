#include <cmath>

template<typename T>
class vector2d
{
public:
    T x,y;
    vector2d(std::initializer_list<T> a) : x{*a.begin()},y{*(a.end()-1)}{}
    vector2d() : x{0},y{0}{}

    vector2d<T> operator+=(vector2d<T> const& b){
        x += b.x;
        y += b.y;
        return *this;
    }

    vector2d<T> operator+(vector2d<T> const& b){
        return vector2d{x+b.x, y+b.y};
    }

    vector2d<T> operator-=(vector2d<T> const& b){
        x -= b.x;
        y -= b.y;
        return *this;
    }

    vector2d<T> operator-(vector2d<T> const& b){
        return vector2d{x-b.x, y-b.y};
    }

    vector2d<T> operator*=(T const& a){
        x = x*a;
        y = y*a;
        return *this;
    }

    template<typename A>
    vector2d<T> operator*(A const& a){
        return vector2d{a*x, a*y};
    }

    vector2d<T> operator/=(T const& a){
        x = x/a;
        y = y/a;
        return *this;
    }

    friend bool operator==(vector2d<T> const& a,vector2d<T> const& b){
        if (a.x == b.x and a.y == b.y){return true;}
        else{return false;}
    }

    friend bool operator!=(vector2d<T> const& a,vector2d<T> const& b){
        if (a.x == b.x and a.y == b.y){return false;}
        else{return true;}
    }

    friend std::ostream& operator<<(std::ostream& o, vector2d<T>const& v){
        o << "[" << v.x << ", " << v.y << "]";
        return o;
    }

    friend std::istream& operator>>(std::istream& i, vector2d<T>& v){
        i >> v.x;
        i >> v.y;
        return i;
    }

    template< typename A>
    friend vector2d<T> operator*(A const& a, vector2d<T> const& v){
        return vector2d{a*v.x, a*v.y};
    }

    friend vector2d<T> normalize(vector2d<T> const& v){
        return vector2d{v.x/length(v), v.y/length(v)};
    }
};

template<typename T>
T dot(vector2d<T> const& a, vector2d<T> const& b){
    return a.x*b.x + a.y*b.y;
}

template<typename T>
T length(vector2d<T> const& v){
    return sqrt(v.x*v.x + v.y*v.y);
}

template<typename T>
T sqlength(vector2d<T> const& v){
    return v.x*v.x + v.y*v.y;
}
