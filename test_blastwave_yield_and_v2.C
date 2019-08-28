//
// Implementation of the blast wave v2 formula
// Klaus Reygers, August 2019
//

// numerator of the blastwave v2 formula
double v2_numerator(const double *x, const double *p) {

    // integration variables
    double rHat = x[0];
    double PhiHat = x[1];

    // parameters
    double pt = p[0];
    double m = p[1];
    double T = p[2];
    double rho0 = p[3];
    double rho2 = p[4];
    double RxOverRy = p[5];

    return rHat *
           TMath::BesselI(
               2,
               (pt * TMath::SinH(rHat * (rho0 + rho2 * TMath::Cos(2 * TMath::ATan(RxOverRy * TMath::Tan(PhiHat)))))) /
                   T) *
           TMath::BesselK(
               1, (TMath::Sqrt(TMath::Power(m, 2) + TMath::Power(pt, 2)) *
                   TMath::CosH(rHat * (rho0 + rho2 * TMath::Cos(2 * TMath::ATan(RxOverRy * TMath::Tan(PhiHat)))))) /
                      T) *
           TMath::Cos(2 * TMath::ATan(RxOverRy * TMath::Tan(PhiHat)));
}

// denominator of the blastwave v2 formula
double v2_denominator(const double *x, const double *p) {

    // integration variables
    double rHat = x[0];
    double PhiHat = x[1];

    // parameters
    double pt = p[0];
    double m = p[1];
    double T = p[2];
    double rho0 = p[3];
    double rho2 = p[4];
    double RxOverRy = p[5];

    return rHat *
           TMath::BesselI(
               0,
               (pt * TMath::SinH(rHat * (rho0 + rho2 * TMath::Cos(2 * TMath::ATan(RxOverRy * TMath::Tan(PhiHat)))))) /
                   T) *
           TMath::BesselK(
               1, (TMath::Sqrt(TMath::Power(m, 2) + TMath::Power(pt, 2)) *
                   TMath::CosH(rHat * (rho0 + rho2 * TMath::Cos(2 * TMath::ATan(RxOverRy * TMath::Tan(PhiHat)))))) /
                      T);
}

void blastwave_yield_and_v2(const double &pt, const double &m, const double &T, const double &rho0, const double &rho2,
                            const double &RxOverRy, double &inv_yield, double &v2) {

    double pars[6] = {pt, m, T, rho0, rho2, RxOverRy};

    // wrapper function for numerical integration of v2 numerator
    ROOT::Math::WrappedParamFunction<> w_v2_num(&v2_numerator, 2, 6);
    w_v2_num.SetParameters(pars);
    ROOT::Math::AdaptiveIntegratorMultiDim ig_num;
    ig_num.SetFunction(w_v2_num);
    ig_num.SetRelTolerance(0.0001);

    // wrapper function for numerical integration of v2 denominator
    ROOT::Math::WrappedParamFunction<> w_v2_den(&v2_denominator, 2, 6);
    w_v2_den.SetParameters(pars);
    ROOT::Math::AdaptiveIntegratorMultiDim ig_den;
    ig_den.SetFunction(w_v2_den);
    ig_den.SetRelTolerance(0.0001);

    // integration range
    double xmin[2] = {0., 0.};
    double xmax[2] = {1., 2. * TMath::Pi()};

    // integrate
    double v2_num = ig_num.Integral(xmin, xmax);
    double v2_den = ig_den.Integral(xmin, xmax);
    if (v2_den != 0) {
        v2 = v2_num / v2_den;
    } else {
        cout << "WARNING: v2 denominator zero!!!" << endl;
    }

    // calculate invariant yield 1/(2 pi pt) dN/(dpt dy) (with arbitrary constant pre-factor)
    inv_yield = TMath::Sqrt(m * m + pt * pt) * v2_den;
}

// main function:
// show usage of blast wave v2 function
void test_blastwave_yield_and_v2() {

    double v2 = 0;
    double inv_yield = 0;

    // blast wave parameters
    const double pt = 1;         // in GeV
    const double m = 0.14;       // in GeV
    const double T = 0.14;       // in GeV
    const double rho0 = 0.925;   // transverse rapidity
    const double rho2 = 0.15;    // azimuthal modulation of transverse rapidity
    const double RxOverRy = 0.9; // ratio of the radii Rx and Ry of the freeze-out ellipse in the transverse plane

    blastwave_yield_and_v2(pt, m, T, rho0, rho2, RxOverRy, inv_yield, v2);

    cout << "v2 = " << v2 << endl;
    cout << "invariant yield = " << inv_yield << endl;
}