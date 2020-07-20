#ifndef STATISTICAL_MODEL
#define STATISTICAL_MODEL

Double_t number_density_analytical(Double_t m, Double_t g, Double_t T, Double_t mu, Double_t sign, Int_t n_terms = 20) {

    // caution: at the moment only for mu = 0!!!

    // [0]: degenerary g
    // [1]: sign, +1 for fermions, -1 for bosons
    // [2]: chemical potential mu
    // [3]: temperature T
    // [4]: mass

    Int_t s = -sign; // fermions have alternating sign in the sum

    Double_t sum_k = 0;

    // n_terms = 1: Boltzmann approximation 
    for (Int_t k = 1; k < n_terms + 1; ++k) {
        sum_k += TMath::Power(s, k + 1) / k * m * m * TMath::BesselK(2, k * m / T);
    }

    return g / (2. * TMath::Pi() * TMath::Pi()) * T * sum_k;
}

Double_t number_density(Double_t m, Double_t spinType, Double_t T, Double_t mu, Int_t n_terms = 20) {

    Double_t sign = 0;
    if (int(spinType) % 2 == 0) {
        //  fermions
        sign = +1;
    } 
    else if (int(spinType) % 2 == 1) {
        // bosons
        sign = -1;
    }

    return number_density_analytical(m, spinType, T, mu, sign, n_terms);

}

#endif


