#include "GF1024.h"
#include <fstream>
#include <cassert>


/********* STATICKE PROMENLJIVE *********/
long GF1024::stepeni[1023]; // niz vrednosti stepena
long GF1024::logaritmi[1024]; // niz vrednosti logaritama


long log(GF1024 n)
{
    return GF1024::logaritmi[n.koef];
}


GF1024::inicijalizator::inicijalizator()
{
    // Funkcija popunjava tabele stepena i logaritama

    long k(1); // trenutni koeficijent
    long const maska(1023); // maska: (dec) 1023 = (bin) 1111111111
    long const generator(1033); // generatorski polinom  z^10 + z^3 + 1
                                // (bin) 10000001001 = (dec) 1033

    stepeni[0] = k;
    logaritmi[k] = 0;
    logaritmi[0] = -1;

    std::cout << "\nPopunjavanje tablica stepena i logaritama...\n";

    for (long i = 1; i < 1023; ++i)
    {
        k <<= 1; // shiftujemo k za jedno mesto u levo
        if (k > maska)
        {
            k ^= generator; // ako je k > mask, radimo XOR sa maskom
        }
        stepeni[i] = k;
        logaritmi[k] = i;
    }
}

GF1024::inicijalizator GF1024::init;


void GF1024::save_tabs()
{
    // Funkcija snima tabele stepena i logaritama u fajlove

    std::fstream f;
    f.open("stepeni1024.txt", std::fstream::out);
    for (int i = 0; i < 1023; ++i)
    {
        f << i << "  " << GF1024::stepeni[i] << std::endl;
    }
    f.close();

    f.open("logaritmi1024.txt", std::fstream::out);
    for (int i = 0; i < 1024; ++i)
    {
        f << i << "  " << GF1024::logaritmi[i] << std::endl;
    }
    f.close();
}

GF1024 GF1024::get_power(long i)
{
    // Funkcija vraca vrednost zeljenog stepena

    if ((i < 0) || (i > 1022))
    {
        i %= 1023;
        if (i < 0)
            i += 1023;
    }
        
    return stepeni[i];
}

GF1024 GF1024::convert_power_to_GF(long power)
{
    // Funkcija sluzi za konvertovanje celog broja u GF1024 koeficijent.
    // 0 ostaje 0, -1 se konvertuje u ceo broj 1 (ne stoji uz alfa),
    // ostali brojevi se citaju iz tablice stepena, jer predstavljaju
    // koeficijente (stepene) polinoma (tj. alfa^stepen).

    if (power == -1) return 0;
    else return GF1024::stepeni[power];
}


GF1024 operator+ (const GF1024& left, const GF1024& right)
{
    // Sabiranje u polju se radi XOR operacijom

    GF1024 res;

    res.koef = left.koef ^ right.koef;
    return res;
}

GF1024 operator- (const GF1024& left, const GF1024& right)
{
    // Oduzimanje je isto kao i sabiranje

    return left + right;
}

GF1024 operator* (const GF1024& left, const GF1024& right)
{
    // Mnozenje se svodi na sabiranje stepena, 
    // tj. logaritmovanih vrednosti koeficijenata

    if ((left.koef == 0) || (right.koef == 0))
    {
        // ako je i jedan koeficijent 0, i rezulat je 0
        return 0;
    }

    long k = log(left) + log(right);

    if (k >= 1023)
    {
        // ako vrednost prelazi 1023, oduzimamo je od rezultata
        k -= 1023;
    }

    return GF1024::get_power(k);  // vracamo dobijenu vrednost u stepene
}

GF1024 operator/ (const GF1024& left, const GF1024& right)
{
    // deljenje se svodi na oduzimanje stepena

    assert(right != 0);

    if (left.koef == 0)
    {
        // ako je deljenik 0, i rezultat je 0
        return 0;
    }

    long k = log(left) - log(right);

    if (k < 0)
    {
        k += 1023;
    }

    return GF1024::get_power(k);
}

bool operator== (const GF1024& left, const GF1024& right)
{
    return left.koef == right.koef;
}

bool operator!= (const GF1024& left, const GF1024& right)
{
    return !(left == right);
}

std::ostream& operator<< (std::ostream& stream, const GF1024& p)
{
    return stream << p.koef;
}

GF1024::operator long() const
{
    return koef;
}

void GF1024::printGF(bool print_zeros) const
{
    if (log(*this) == -1)
    {
        if (print_zeros == true)
            std::cout << "0";
    }
    else if (log(*this) == 0)
        std::cout << "1";
    else if (log(*this) == 1)
        std::cout << "a";
    else
        std::cout << "a^" << log(*this);
}
