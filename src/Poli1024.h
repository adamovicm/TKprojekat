#ifndef POLI1024_H
#define POLI1024_H

#include "GF1024.h"
#include <vector>


using std::vector;

class Poli1024 : public vector<GF1024>
{

public:
    Poli1024(GF1024 k = 0) : vector<GF1024> (k == 0 ? 0 : 1, k) {}
    Poli1024(long nr, GF1024 value) : vector<GF1024>(nr, value) {}
    Poli1024(const vector<GF1024>& v) : vector<GF1024>(v) {}
    
    long power() const;
    void normalize();

private:
    static void divide_poli(Poli1024& u, Poli1024& v);

public:
    friend Poli1024 operator+ (const Poli1024& left, const Poli1024& right);
    friend Poli1024 operator* (const Poli1024& left, const Poli1024& right);
    friend Poli1024 operator/ (Poli1024 left, Poli1024 right);
    friend Poli1024 operator% (Poli1024 left, Poli1024 right);
    friend bool operator== (const Poli1024& left, const Poli1024& right);
    friend bool operator!= (const Poli1024& left, const Poli1024& right);

    GF1024 horner(const GF1024& Z);

    void print_poli() const;

};

#endif
