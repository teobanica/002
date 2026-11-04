#pragma once
#include <numeric>
#include <cassert>
namespace boost {
template<class T> class rational {
 T n,d; void norm(){assert(d);if(d<0){n=-n;d=-d;}T g=std::gcd(n,d);n/=g;d/=g;}
 public:rational(T a=0,T b=1):n(a),d(b){norm();} T numerator()const{return n;} T denominator()const{return d;}
 rational operator-()const{return rational(-n,d);}
 rational& operator+=(const rational&b){n=n*b.d+b.n*d;d*=b.d;norm();return *this;}
 rational& operator-=(const rational&b){return *this+=-b;}
 rational& operator*=(const rational&b){n*=b.n;d*=b.d;norm();return *this;}
 rational& operator/=(const rational&b){n*=b.d;d*=b.n;norm();return *this;}
 friend rational operator+(rational a,const rational&b){return a+=b;}
 friend rational operator-(rational a,const rational&b){return a-=b;}
 friend rational operator*(rational a,const rational&b){return a*=b;}
 friend rational operator/(rational a,const rational&b){return a/=b;}
 friend bool operator==(const rational&a,const rational&b){return a.n==b.n&&a.d==b.d;}
 friend bool operator!=(const rational&a,const rational&b){return !(a==b);}
 friend bool operator<(const rational&a,const rational&b){return (__int128)a.n*b.d<(__int128)b.n*a.d;}
};}
