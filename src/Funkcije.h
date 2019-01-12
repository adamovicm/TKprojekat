#ifndef FUNKCIJE_H
#define FUNKCIJE_H

#include "Constants.h"
#include "GF1024.h"
#include "Poli1024.h"

using std::vector;

Poli1024 Generator();
Poli1024 Coder(const Poli1024& inf_seq, const Poli1024& gen);

std::vector<GF1024> FFTdef33(vector<GF1024>& niz);

std::vector<GF1024> GDFFT3x11(const vector<GF1024>& niz);
Poli1024 GTFFT(const Poli1024& niz);
std::vector<GF1024> Sindrom(const Poli1024& received_seq);
Poli1024 BerlekampMassey(const vector<GF1024>& sindrom);
vector<long> Find_zeros(Poli1024 Sigma);
Poli1024 Forney(const Poli1024& Sigma, const vector<GF1024>& Sindrom, const vector<long>& zero_places);

Poli1024 FFTdef(Poli1024& niz);
std::vector<GF1024> SindromDef(Poli1024& primljena_sekvenca, long tD, long N);


#endif
