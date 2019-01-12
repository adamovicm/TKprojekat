#ifndef GF1024_H
#define GF1024_H

#include <iostream>

class GF1024
{
private:
    long koef;   // koeficijenti
    static long stepeni[1023];  // tabela stepena
    static long logaritmi[1024];  // tabela logaritama

    class inicijalizator
    {
    public:
        inicijalizator();
    };

    static inicijalizator init;

public:
    GF1024(long k = 0) : koef(k & 1023) {}

    static void save_tabs();
    static GF1024 get_power(long i);
    static GF1024 convert_power_to_GF(long power);

    friend long log(GF1024 n);

    friend GF1024 operator+ (const GF1024& left, const GF1024& right);
    friend GF1024 operator- (const GF1024& left, const GF1024& right);
    friend GF1024 operator* (const GF1024& left, const GF1024& right);
    friend GF1024 operator/ (const GF1024& left, const GF1024& right);

    friend bool operator== (const GF1024& left, const GF1024& right);
    friend bool operator!= (const GF1024& left, const GF1024& right);

    friend std::ostream& operator<< (std::ostream& stream, const GF1024& p);

    explicit operator long() const;

    void printGF(bool print_zeros = false) const;
};

#endif // !GF1024_H
