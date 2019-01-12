#include <iostream>
#include <vector>
#include "GF1024.h"
#include "Poli1024.h"
#include "Funkcije.h"
#include "Constants.h"

using std::cout;
using std::endl;
using std::vector;

int main()
{
    // PO POTREBI: snimanje tabela stepena i logaritama u fajlove
    //GF1024::snimiTabeleUFajl(); // tabele se snimaju u *.txt fajlove

    // Generator
    cout << "\n\nIzracunavanje generatorskog polinoma:\n";
    Poli1024 gen = Generator();  // generatorski polinom
    cout << "\ng(Z) = ";
    gen.print_poli();
    cout << endl;

    /*************** INFORMACIONANA SEKVENCA ***************/
    // zadati inf sekvencu
        
    Poli1024 inf_seq(K, 2);  // informaciona sekvenca

    cout << "\n\nKodovanje...\n";
    Poli1024 code_word = Coder(inf_seq, gen);  // kodna rec
    Poli1024 received_seq = code_word;  // kodna rec nakon dodatih gresaka

    /*************** DODAVANJE GRESAKA ***************/
    // dodati greske polinomu

    received_seq[100] = 5;
    received_seq[202] = 5;
    received_seq[333] = 5;
    received_seq[401] = 5;
    received_seq[550] = 5;
    //received_seq[50] = 5;
    //received_seq[55] = 5;

    /*********** TESTOVI **********/
    //Poli1024 ftDef, ftGT, div;
    //div = dividePoli(kodnaRec, gen, 'r');
    //ftDef = FFTdef(primljena_sekvenca);
    //ftGT = GTFFT(primljena_sekvenca);
    //if (ftDef == ftGT) std::cout << "\nok!\n";
    //std::vector<GF1024> sinDef(2 * tD, 0);
    //sinDef = SindromDef(received_seq, tD, N);


    /*********** PRORACUN **********/

    // Sindrom
    cout << "\n\nIzracunavanje sindroma:\n\n";
    vector<GF1024> sindrom = Sindrom(received_seq);  // sindrom
    for (long i = 0; i < 2 * tD; ++i)
    {
        cout << "S" << i + 1 << " = ";
        sindrom[i].printGF();
        cout << endl;
    }
    
    // TEST
    //if (sinDef == sindrom) cout << "\n\nJEDNAKI SU!\n\n";

    // Berlekamp-Massey
    cout << "\n\nBerlekamp Massey-ev algoritam:\n";
    Poli1024 Sigma = BerlekampMassey(sindrom);  // polinom lokacije gresaka
    cout << "\nSigma(Z) = ";
    Sigma.print_poli();
    cout << endl;
    // Stepen polinoma Sigma ne sme biti veci od tD
    // (ne moze biti vise gresaka nego sto algoritam moze da ispravi)
    if (Sigma.power() > tD)
    {
        cout << "\n\nStepen polinoma Sigma > tD!\n\n";
        cout << "Greske se ne mogu ispraviti.\n\n";
        cout << "\n\nPritisnite [enter] za kraj\n";
        std::cin.ignore();
        return 0;
    }


    // Trazenje nula
    cout << "\n\nFaktorisanje polinoma Sigma:\n";
    vector<long> zero_places = Find_zeros(Sigma);  // mesta gresaka
    cout << "\nBroj pronadjenih gresaka: " << zero_places.size() << "\n";
    cout << "Mesta gresaka: ";
    for (long i = 0; i < zero_places.size(); ++i)
    {
        cout << N - zero_places[i] << " ";
    }
    cout << endl;
    // Broj nula polinoma Sigma i njegov stepen moraju biti isti!
    if (zero_places.size() != Sigma.power())
    {
        cout << "\n\nStepen polinoma Sigma nije jednak broju njegovih nula!\n\n";
        cout << "Greske se ne mogu ispraviti.\n\n";
        cout << "\n\nPritisnite [enter] za kraj\n";
        std::cin.ignore();
        return 0;
    }
        

    // Forney-ev algoritam
    cout << "\n\nForney-ev algoritam:\n";
    Poli1024 errors = Forney(Sigma, sindrom, zero_places);  // vrednosti gresaka
    cout << "\ne(Z) = ";
    errors.print_poli();

    // Ispravka gresaka
    Poli1024 x = received_seq + errors;  // ispravljen polinom

    // TEST
    //x[2] = 0;

    if (x == code_word)
        cout << "\n\n\nSekvenca je uspesno ispravljena!\n\n";
    else
        cout << "\n\n\nSekvenca NIJE ispravljena!!\n\n";

    cout << "\n\nPritisnite [enter] za kraj\n";
    std::cin.ignore();
    
    return 0;
}