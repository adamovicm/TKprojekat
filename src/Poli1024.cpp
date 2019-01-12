#include "Poli1024.h"
#include <algorithm>
#include <cassert>

long Poli1024::power() const
{
    // Funkcija trazi stepen polinoma
    if (empty())
        return -1;

    for (unsigned long i = size() - 1; i > 0; --i)
    {
        if ((*this)[i] != 0)
        {
            return i;
        }
    }

    return 0;
}


void Poli1024::normalize()
{
    // Funkcija normalizuje polinom, tj. uklanja
    // nule sa njegovog kraja

    while (!empty() && (back() == 0))
    {
        pop_back();
    }
}


void Poli1024::divide_poli(Poli1024& u, Poli1024& v)
{
    // Funkcija deli dva polinoma i nalazi dva rezultata:
    //  - ostatak koji ce na kraju funkcije biti uneti parametar u
    //  - kolicnik koji ce biti vracen kao uneti parametar v

    // ako je delilac prazan polinom ili polinom nultog stepena
    // a da mu je slobodan clan = 0, deljenje nije definisano
    assert(!v.empty() && !(v.power() == 0 && v[0] == 0));

    if (u.empty() || (u.power() == 0 && u[0] == 0))
    {
        v.clear();
        return;
    }

    // stepeni polinoma u i v
    long upower = u.power();
    long vpower = v.power();

    Poli1024 q(upower - vpower + 1, 0);

    while (upower >= vpower)
    {
        q[upower - vpower] = u[upower] / v[vpower];

        // mnozenje dobijenog kolicnika sa deliocem
        Poli1024 tmp(upower + 1, 0);
        for (long i = vpower; i >= 0; --i)
        {
            tmp[i + upower - vpower] = v[i] * q[upower - vpower];
        }

        // racunanje ostatka sto je zapravo i novi deljenik
        for (long i = upower; i >= upower - vpower; --i)
        {
            u[i] = u[i] + tmp[i];
        }

        u.normalize();
        upower = u.power();
    }

    v = q;
}


Poli1024 operator+ (const Poli1024& left, const Poli1024& right)
{
    Poli1024 res;

    long lsize = left.size();
    long rsize = right.size();

    // U slucaju da su polinomi razlicitih duzina, sabiramo
    // njihove elemente do kraja polinoma manjeg stepena
    for (long i = 0; i < std::min(lsize, rsize); ++i)
    {
        res.push_back(left[i] + right[i]);
    }

    // Ako su polinomi razlicitih duzina, u rezultat dodajemo ostatak
    // polinoma veceg stepena
    if (lsize < rsize)
    {
        for (long i = lsize; i < rsize; ++i)
        {
            res.push_back(right[i]);
        }
    }
    else if (rsize < lsize)
    {
        for (long i = rsize; i < lsize; ++i)
        {
            res.push_back(left[i]);
        }
    }

    res.normalize();
    
    return res;
}

Poli1024 operator* (const Poli1024& left, const Poli1024& right)
{
    // rezultat ce biti maksimalne duzine kao i zbir stepena kolicnika
    if (left.empty() || right.empty())
        return Poli1024();

    long lpower = left.power();
    long rpower = right.power();

    Poli1024 res((lpower + rpower + 1), 0);

    for (long l = 0; l <= lpower; ++l)
    {
        for (long r = 0; r <= rpower; ++r)
        {
            long stepen = l + r;
            res[stepen] = res[stepen] + (left[l] * right[r]);
        }
    }

    res.normalize();
    
    return res;
}

Poli1024 operator/ (Poli1024 left, Poli1024 right)
{
    Poli1024::divide_poli(left, right);
    return right;
}

Poli1024 operator% (Poli1024 left, Poli1024 right)
{
    Poli1024::divide_poli(left, right);
    return left;
}

bool operator==(const Poli1024& left, const Poli1024& right)
{
    //assert(!right.empty());
    //assert(!(right.power() == 0 && right[0] == 0));

    long lpower = left.power();
    long rpower = right.power();

    // ako su polinomi razlicitih stepena => nisu jednaki
    if (lpower != rpower)
        return false;

    for (long i = 0; i <= lpower; ++i)
    {
        if (left[i] != right[i])
            return false;
    }
    
    return true;
}

bool operator!=(const Poli1024& left, const Poli1024& right)
{
    return !(left == right);
}

GF1024 Poli1024::horner(const GF1024& Z)
{
    // Funkcija trazi vrednost polinoma u tacki Z
    // koristeci Hornerovu semu

    if (Z == 0)
        return 0;

    GF1024 result(back());
    for (long n = power() - 1; n >= 0; --n)
    {
        result = result * Z + (*this)[n];
    }

    return result;
}


void Poli1024::print_poli() const
{
    // Funkcija ispisuje polinom (stvarne stepene koeficijenata)

    // prvi - da li je prvi nenulti element u pitanju, ispred njega ne dodajemo +
    bool prvi(true);
    
    for (unsigned long i = 0; i < size(); ++i)
    {
        if ((*this)[i] != 0)
        {
            // ako nije prvi nenulti element, prvo dodajemo +
            if (prvi != true)
                std::cout << " + ";

            // ako je nulti element, ispisujemo samo njega
            if (i == 0)
                (*this)[i].printGF();

            // ako je prvi element, dodajemo Z nakon njega
            else if (i == 1)
            {
                // ako je element ==1, ispisujemo samo Z
                if ((*this)[i] == 1)
                    std::cout << "Z";
                else
                {
                    (*this)[i].printGF();
                    std::cout << ".Z";
                }
            }

            // ako nije ni nulti ni prvi element, ispisujemo ga i dodajemo Z^i
            else
            {
                if ((*this)[i] == 1)
                    std::cout << "Z^" << i;
                else
                {
                    (*this)[i].printGF();
                    std::cout << ".Z^" << i;
                }
            }

            prvi = false;
        }
    }
}