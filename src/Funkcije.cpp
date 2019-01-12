#include "Funkcije.h"
#include <vector>
#include <fstream>
#include <cassert>

using std::vector;

Poli1024 Generator()
{
    // Funkcija pravi generatorski niz
    // g(z) = (z + alfa)(z + alfa^2)(z + alfa^3)...

    // g se postavlja na (z + alfa)
    Poli1024 g(2,0);
    g[0] = 2;
    g[1] = 1;

    // tmp na pocetku postavljamo na (z)
    // a u petlji dodajemo alfa^2, alfa^3...
    Poli1024 tmp(2, 0);
    tmp[1] = 1;

    for (long i = 2; i <= 2 * tD; ++i)
    {
        tmp[0] = GF1024::get_power(i);
        g = g * tmp;
    }

    return g;
}

Poli1024 Coder(const Poli1024& inf_seq, const Poli1024& gen)
{
    // Funkcija koduje informativnu sekvencu i vraca kodnu rec

    Poli1024 tmp = inf_seq;
    Poli1024::iterator iter = tmp.begin();


    tmp.insert(iter, N - K, GF1024(0));

    Poli1024 code_word = tmp % gen;
    return code_word + tmp;
}


std::vector<GF1024> FFTdef33(std::vector<GF1024>& niz)
{
    // Furijeova transofrmacija "po definiciji" za niz duzine 33.
    // Za svrhe testiranja ispravnosti rezultata

    long N(1023), n(33);
    std::vector<GF1024> rezultat(33,0);

    for (long k = 0; k < n; ++k)
    {
        GF1024 Bk;
        for (long i = 0; i < n; ++i)
        {
            GF1024 omega = GF1024::get_power((i * k * N/n) % N);
            Bk = Bk + omega * niz[i];
        }
        rezultat[k] = Bk;
    }
    return rezultat;
}


std::vector<GF1024> GDFFT3x11(std::vector<GF1024>& niz)
{
    // Funkcija trazi Furijeovu transformaciju za niz duzine 33 (3x11)
    // koristeci Good-Thomas algoritam

    long N1N2(33), N1(3), N2(11);

    // 1D -> 2D
    // matrix - matrica za smestanje elemenata niza (1D niz u 2D matricu) dimenzija 3x11
    vector<vector<GF1024>> matrix(N1, vector<GF1024>(N2, 0));
    for (long i = 0; i < N1N2; ++i)
        matrix[i % N1][i % N2] = niz[i];

    // temp_matrix - matrica za privremeno smestanje rezultata FT
    std::vector<std::vector<GF1024>> temp_matrix(N1, std::vector<GF1024>(N2, 0));
    GF1024 Bk, omega;

    // prvo radimo FT po vrstama
    for (long vrsta = 0; vrsta < N1; ++vrsta)
    {
        for (long k = 0; k < N2; ++k)
        {
            Bk = 0;
            for (long i = 0; i < N2; ++i)
            {
                omega = GF1024::get_power(((N / N2) * i * k) % N);
                Bk = Bk + omega * matrix[vrsta][i];
            }
            temp_matrix[vrsta][k] = Bk;
        }
    }

    matrix = temp_matrix;

    // zatim radimo FT po kolonama
    for (long kolona = 0; kolona < N2; ++kolona)
    {
        for (long k = 0; k < N1; ++k)
        {
            Bk = 0;
            for (long i = 0; i < N1; ++i)
            {
                omega = GF1024::get_power(((N / N1) * i * k) % N);
                Bk = Bk + omega * matrix[i][kolona];
            }
            temp_matrix[k][kolona] = Bk;
        }
    }

    // 2D -> 1D
    std::vector<GF1024> result(N1N2, 0);
    for (long k1 = 0; k1 < N1; ++k1)
        for (long k2 = 0; k2 < N2; ++k2)
            result[(N2 * k1 + N1 * k2) % N1N2] = temp_matrix[k1][k2];

    return result;
}


Poli1024 GTFFT(const Poli1024& niz)
{
    // Funkcija trazi Furijeovu transformaciju za niz duzine 1023
    // koristeci Good-Thomas algoritam

    long N1(31), N2(33);

    // 1D -> 2D
    // matrix - matrica za smestanje elemenata niza (1D niz u 2D matricu) dimenzija 3x11
    vector<vector<GF1024>> matrix(N1, vector<GF1024>(N2, 0));
    for (long i = 0; i < N; ++i)
        matrix[i % N1][i % N2] = niz[i];

    // pod transformacije
    for (long i = 0; i < N1; ++i)
        matrix[i] = GDFFT3x11(matrix[i]);
    // ovime su FT po vrstama vec odradjene


    // radimo FT po kolonama
    // temp_matrix - matrica za privremeno smestanje rezultata FT
    vector<vector<GF1024>> temp_matrix(N1, vector<GF1024>(N2, 0));
    GF1024 Bk, omega;
    for (long kolona = 0; kolona < N2; ++kolona)
    {
        for (long k = 0; k < N1; ++k)
        {
            Bk = 0;
            for (long i = 0; i < N1; ++i)
            {
                omega = GF1024::get_power(((N / N1) * i * k) % N);
                Bk = Bk + omega * matrix[i][kolona];
            }
            temp_matrix[k][kolona] = Bk;
        }
    }

    // 2D -> 1D
    Poli1024 result(N, 0);
    for (long k1 = 0; k1 < N1; ++k1)
        for (long k2 = 0; k2 < N2; ++k2)
            result[(N2 * k1 + N1 * k2) % N] = temp_matrix[k1][k2];

    return result;
}

vector<GF1024> Sindrom(const Poli1024& received_seq)
{
    // Funkcija trazi sindrome

    vector<GF1024> sindrom(2 * tD, 0);
    Poli1024 ft = GTFFT(received_seq);

    // Sindrom se dobije iz FT primljene sekvence iz koje uzimamo prvih
    // 2*tD rezultata, prekakajuci nulti.
    for (long i = 0; i < 2 * tD; ++i)
        sindrom[i] = ft[i + 1];        

    return sindrom;
}

Poli1024 BerlekampMassey(const vector<GF1024>& sindrom)
{
    long l(1), m(0), Ll_1(0);
    GF1024 Delta_m(1), Delta_l;

    // inicijalizacija
    Poli1024 Sigma_l_1(1, 1);
    Poli1024 Sigma_m_1(1, 1);
    Poli1024 Sigma_l(1, 1);

    while (l <= 2 * tD)
    {
        // Delta_l
        Delta_l = sindrom[l - 1];
        for (long k = 1; k <= Ll_1; ++k)
            Delta_l = Delta_l + (Sigma_l_1[k] * sindrom[(l - 1) - k]);

        if (Delta_l != 0)
        {
            Poli1024 Sigma_m_1_tmp = Sigma_m_1 * (Delta_l / Delta_m);
            Poli1024::iterator iter = Sigma_m_1_tmp.begin();
            Sigma_m_1_tmp.insert(iter, l - m, 0);
            Sigma_l = Sigma_l_1 + Sigma_m_1_tmp;

            // promena indeksa, formalnog stepena i radnih promenljivih
            // kada je 2 * Ll-1 < l
            if (2 * Ll_1 < l)
            {
                m = l;
                Ll_1 = l - Ll_1;
                Delta_m = Delta_l;
                Sigma_m_1 = Sigma_l_1;
            }

            Sigma_l_1 = Sigma_l;
        }

        ++l;
    }
 
    return Sigma_l;
}

vector<long> Find_zeros(Poli1024 Sigma)
{
    // Algoritam za trazenje nula polinoma Sigma

    // e_place - lokacije gresaka 
    // inv_e_place - logaritam reciprocne vrednosti lokacija gresaka

    // Dopunjavamo polinom Sigma nulama do duzine N
    Sigma.resize(N);

    // Beta - faktorizovani polinom Sigma, za trazenje nula
    Poli1024 Beta = GTFFT(Sigma);

    vector<long> zero_places;
    for (unsigned long i = 0; i < Beta.size(); ++i)
        if (Beta[i] == 0)
            zero_places.push_back(i);

    return zero_places;
}


Poli1024 Forney(const Poli1024& Sigma, const vector<GF1024>& Sindrom, const vector<long>& zero_places)
{
    // Fornijev algoritam za ispravljanje gresaka

    // gama - vrednost greske koja se trenutno racuna
    
    // omega = R_z^(ni+1)  (S(Z)*Sigma(Z)) - evaluator gresaka
    // e - greska

    // SZ - S(Z). Polinom sa sindromima kao koeficijentima. Nema nulti clan
    Poli1024 SZ(2 * tD + 1, 0);
    for (long i = 0; i < 2 * tD; ++i)
        SZ[i + 1] = Sindrom[i];

    // ni - broj pronadjenih nula
    long ni = zero_places.size();

    Poli1024 omega = SZ * Sigma;
    omega.resize(ni + 1);

    // Polinom e ce biti stepena reciprocne vrednosti nule na najnizem mestu
    Poli1024 e(N - zero_places[0] + 1, 0);

    // Petlja prolazi kroz sve greske i smesta vrednosti u e
    for (long k = 0; k < ni; ++k)
    {
        // Bk_1 mesto nule
        GF1024 Bk_1 = GF1024::get_power(zero_places[k]);
        //GF1024 Bk_1 = log(GF1024(zero_places[k]));

        // delilac
        GF1024 denominator(1);
        for (long l = 0; l < ni; ++l)
        {  
            if (l != k)
            {
                // Bl - reciprocna vrednost mesta nule
                GF1024 Bl = GF1024::get_power(N - zero_places[l]);
                //GF1024 Bl = log(GF1024(N - zero_places[k]));
                denominator = denominator * (1 + Bl * Bk_1);
            }   
        }
        assert(denominator != 0);

        // deljenik
        GF1024 numerator(omega.horner(Bk_1));

        e[N - zero_places[k]] = numerator / denominator;
    }

    e.normalize();

    return e;
}



/************ FUNKCIJE "PO DEFINICIJI" **************/

Poli1024 FFTdef(Poli1024& niz)
{
    long N(1023);
    GF1024 omega, Bk;
    Poli1024 rezultat;

    rezultat.assign(N, 0);

    for (long k = 0; k < N; ++k)
    {
        Bk = 0;
        for (long i = 0; i < N; ++i)
        {
            omega = GF1024::get_power((i * k) % N);
            Bk = Bk + omega * niz[i];
        }
        rezultat[k] = Bk;
    }
    return rezultat;
}

std::vector<GF1024> SindromDef(Poli1024& primljena_sekvenca, long tD, long N)
{
    std::vector<GF1024> sindrom(2 * tD, 0);
    GF1024 tmp;

    for (long i = 1; i <= 2 * tD; ++i)
    {
        for (unsigned long j = 0; j < primljena_sekvenca.size(); ++j)
        {
            tmp = GF1024::get_power((i * j) % N);
            sindrom[i - 1] = sindrom[i - 1] + (primljena_sekvenca[j] * tmp);
        }
    }
    return sindrom;
}