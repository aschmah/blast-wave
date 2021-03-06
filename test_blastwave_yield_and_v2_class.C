//
// Implementation of the blast wave v2 formula for an elliptical freeze-out surface
// Klaus Reygers, November 2019
//

class blastwave_yield_and_v2 {

    // function for v2 numerator and denominator
    // fos = freeze-out surface
    // fos1: freeze-out time tf in the lab system: tf = sqrt(tau^2 + z^2)
    // boost: Boltzmann distribution in the cell frame boosted to lab frame (no explicity freeze-out hyper-surface)
    static double v2_fos1_numerator(const double *x, const double *p);
    static double v2_fos1_denominator(const double *x, const double *p);
    static double v2_boost_numerator(const double *x, const double *p);
    static double v2_boost_denominator(const double *x, const double *p);

    // wrapper functions for v2 numerator and denominator
    ROOT::Math::WrappedParamFunction<> w_v2_fos1_num;
    ROOT::Math::WrappedParamFunction<> w_v2_fos1_den;
    ROOT::Math::WrappedParamFunction<> w_v2_boost_num;
    ROOT::Math::WrappedParamFunction<> w_v2_boost_den;

    ROOT::Math::AdaptiveIntegratorMultiDim ig;

  public:
    blastwave_yield_and_v2()
        : w_v2_fos1_num(&blastwave_yield_and_v2::v2_fos1_numerator, 2, 7),
          w_v2_fos1_den(&blastwave_yield_and_v2::v2_fos1_denominator, 2, 7),
          w_v2_boost_num(&blastwave_yield_and_v2::v2_boost_numerator, 2, 7),
          w_v2_boost_den(&blastwave_yield_and_v2::v2_boost_denominator, 2, 7) {}

    void calc_blastwave_yield_and_v2_fos1(const double &pt, const double &m, const double &T, const double &rho0,
                                          const double &rho2, const double &RxOverRy, double &inv_yield, double &v2);

    void calc_blastwave_yield_and_v2_boost(const double &pt, const double &m, const double &T, const double &rho0,
                                           const double &rho2, const double &RxOverRy, double &inv_yield, double &v2);

    ClassDef(blastwave_yield_and_v2, 1)
};

// numerator of the blastwave v2 formula
double blastwave_yield_and_v2::v2_fos1_numerator(const double *x, const double *p) {

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
    double A = p[6]; // arbitrary factor, can be adjusted to improve numerical stability

    double PhiB = TMath::ATan(RxOverRy * TMath::Tan(PhiHat)); // boost angle
    double mt = TMath::Sqrt(m * m + pt * pt);
    double rho = rHat * (rho0 + rho2 * TMath::Cos(2 * PhiB)); // transverse rapidity
    double xip = pt * TMath::SinH(rho) / T;
    double xim = mt * TMath::CosH(rho) / T;

    return A * rHat * TMath::BesselI(2, xip) * TMath::BesselK(1, xim) * TMath::Cos(2 * PhiB);
}

// denominator of the blastwave v2 formula
double blastwave_yield_and_v2::v2_fos1_denominator(const double *x, const double *p) {

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
    double A = p[6]; // arbitrary factor, can be adjusted to improve numerical stability

    double PhiB = TMath::ATan(RxOverRy * TMath::Tan(PhiHat)); // boost angle
    double mt = TMath::Sqrt(m * m + pt * pt);
    double rho = rHat * (rho0 + rho2 * TMath::Cos(2 * PhiB)); // transverse rapidity
    double xip = pt * TMath::SinH(rho) / T;
    double xim = mt * TMath::CosH(rho) / T;

    return A * rHat * TMath::BesselI(0, xip) * TMath::BesselK(1, xim);
}

void blastwave_yield_and_v2::calc_blastwave_yield_and_v2_fos1(const double &pt, const double &m, const double &T,
                                                              const double &rho0, const double &rho2,
                                                              const double &RxOverRy, double &inv_yield, double &v2) {

    // blast wave parameters:
    // the last number (par[6]) is an arbitrary normalization which we will adjust
    // in order to have a good numerical stability for different masses and pt
    double pars[7] = {pt, m, T, rho0, rho2, RxOverRy, 1.};

    // determine scale factor which ensures good numerical stability
    const double xsf[2] = {1., 0.};
    double sf = v2_fos1_denominator(xsf, pars);
    pars[6] = 1. / sf;

    // set parameters
    w_v2_fos1_num.SetParameters(pars);
    w_v2_fos1_den.SetParameters(pars);

    // define integrator
    // ROOT::Math::AdaptiveIntegratorMultiDim ig;
    ig.SetRelTolerance(1e-6);

    // integration range
    double xmin[2] = {0., 0.};
    double xmax[2] = {1., 2. * TMath::Pi()};

    // integrate
    ig.SetFunction(w_v2_fos1_num);
    double v2_num = ig.Integral(xmin, xmax);
    // if (ig_num.Status() != 0) cout << ig_num.Status() << endl;

    ig.SetFunction(w_v2_fos1_den);
    double v2_den = ig.Integral(xmin, xmax);
    // if (ig_den.Status() != 0) cout << ig_den.Status() << endl;

    // cout << pt << " " << v2_den << endl;
    // cout << pt << " " << inv_yield << endl;

    if (v2_den != 0) {
        v2 = v2_num / v2_den;
    } else {
        cout << "WARNING: v2 denominator zero!!!" << endl;
    }

    inv_yield = sf * TMath::Sqrt(m * m + pt * pt) * v2_den;
}

// numerator of the blastwave v2 formula
double blastwave_yield_and_v2::v2_boost_numerator(const double *x, const double *p) {

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
    double A = p[6]; // arbitrary factor, can be adjusted to improve numerical stability

    double PhiB = TMath::Pi() * TMath::Floor((PhiHat + TMath::Pi() / 2.) / TMath::Pi()) +
                  TMath::ATan(RxOverRy * TMath::Tan(PhiHat)); // boost angle
    double rho = rHat * (rho0 + rho2 * TMath::Cos(2 * PhiB)); // transverse rapidity
    double mt = TMath::Sqrt(m * m + pt * pt);
    double xip = pt * TMath::SinH(rho) / T;
    double xim = mt * TMath::CosH(rho) / T;

    return A * rHat * TMath::Cos(2 * PhiB) *
           (TMath::BesselI(2, xip) * (2 * T * TMath::BesselK(0, xim) + mt * TMath::BesselK(1, xim) * TMath::CosH(rho)) -
            pt * TMath::BesselI(1, xip) * TMath::BesselK(0, xim) * TMath::SinH(rho));
}

// denominator of the blastwave v2 formula for freeze-out surface (fos) 2
double blastwave_yield_and_v2::v2_boost_denominator(const double *x, const double *p) {

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
    double A = p[6]; // arbitrary factor, can be adjusted to improve numerical stability
    double PhiB = TMath::Pi() * TMath::Floor((PhiHat + TMath::Pi() / 2.) / TMath::Pi()) +
                  TMath::ATan(RxOverRy * TMath::Tan(PhiHat)); // boost angle
    double rho = rHat * (rho0 + rho2 * TMath::Cos(2 * PhiB)); // transverse rapidity
    double mt = TMath::Sqrt(m * m + pt * pt);
    double xip = pt * TMath::SinH(rho) / T;
    double xim = mt * TMath::CosH(rho) / T;

    return A * rHat *
           (mt * TMath::BesselI(0, xip) * TMath::BesselK(1, xim) * TMath::CosH(rho) -
            pt * TMath::BesselI(1, xip) * TMath::BesselK(0, xim) * TMath::SinH(rho));
}

void blastwave_yield_and_v2::calc_blastwave_yield_and_v2_boost(const double &pt, const double &m, const double &T,
                                                               const double &rho0, const double &rho2,
                                                               const double &RxOverRy, double &inv_yield, double &v2) {

    // blast wave parameters:
    // the last number (par[6]) is an arbitrary normalization which we will adjust
    // in order to have a good numerical stability for different masses and pt
    double pars[7] = {pt, m, T, rho0, rho2, RxOverRy, 1.};

    // determine scale factor which ensures good numerical stability
    const double xsf[2] = {1., 0.};
    double sf = v2_boost_denominator(xsf, pars);
    // cout << "debugging: sf = " << sf << endl;
    // double tmp = v2_boost_numerator(xsf, pars);
    // cout << "debugging: tmp = " << tmp << endl;
    pars[6] = 1. / sf;

    // set parameters
    w_v2_boost_num.SetParameters(pars);
    w_v2_boost_den.SetParameters(pars);

    // define integrator
    // ROOT::Math::AdaptiveIntegratorMultiDim ig;
    ig.SetRelTolerance(1e-6);

    // integration range
    double xmin[2] = {0., 0.};
    double xmax[2] = {1., 2. * TMath::Pi()};

    // integrate
    ig.SetFunction(w_v2_boost_num);
    double v2_num = ig.Integral(xmin, xmax);

    ig.SetFunction(w_v2_boost_den);
    double v2_den = ig.Integral(xmin, xmax);
    // cout << "v2_den boost: " << v2_den << endl;

    // // cout << pt << " " << v2_den << endl;
    // // cout << pt << " " << inv_yield << endl;

    if (v2_den != 0) {
        v2 = v2_num / v2_den;
    } else {
        cout << "WARNING: v2 denominator zero!!!" << endl;
    }

    inv_yield = sf * v2_den;
}

// main function:
// show usage of blast wave v2 function
void test_blastwave_yield_and_v2_class() {

    blastwave_yield_and_v2 bw;

    //
    // test 1: single pt point
    //
    double v2 = 0;
    double inv_yield = 0;

    // blast wave parameters
    const double pt = 1;         // in GeV
    const double m = 0.14;       // in GeV
    const double T = 0.14;       // in GeV
    const double rho0 = 0.925;   // transverse rapidity
    const double rho2 = 0.15;    // azimuthal modulation of transverse rapidity
    const double RxOverRy = 0.9; // ratio of the radii Rx and Ry of the freeze-out ellipse in the transverse plane

    bw.calc_blastwave_yield_and_v2_fos1(pt, m, T, rho0, rho2, RxOverRy, inv_yield, v2);

    cout << "Test 1: inv. yield and v2 at a fixed point (for comparison with Mathematica result)" << endl;
    cout << "v2 = " << v2 << endl;
    cout << "invariant yield = " << inv_yield << endl;

    //
    // test 2: plot v2 vs pt
    //
    const double T_BW = 0.102846;
    const double Rho0_BW = 1.14149;
    const double Rho2_BW = 0.0652361;
    const double RxOverRy_BW = 0.824434;

    // const double m_BW = 9.460;
    // const double m_BW = 3.096;
    // const double m_BW = 0.14;
    const double m_BW = 0.938;

    const double pt_min = 0;
    const double pt_max = 16.;
    const double dpt = 0.1;

    const int i_max = (pt_max - pt_min) / dpt + 1;
    TGraph giy(i_max);
    TGraph gv2(i_max);

    cout << endl << "Test 2: Performance test (plot v2 vs pt)" << endl;

    TStopwatch t;
    t.Start();

    for (int i = 0; i < i_max; ++i) {
        double pt_BW = pt_min + i * dpt + 0.2 * dpt;
        double inv_yield_BW = 0;
        double v2_BW = 0;
        bw.calc_blastwave_yield_and_v2_fos1(pt_BW, m_BW, T_BW, Rho0_BW, Rho2_BW, RxOverRy_BW, inv_yield_BW, v2_BW);
        giy.SetPoint(i, pt_BW, inv_yield_BW);
        gv2.SetPoint(i, pt_BW, v2_BW);
    }

    cout << "End" << endl;

    t.Stop();
    t.Print();

    // plot v2 vs pt
    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1");
    TH2F *fr = new TH2F("frame", "blast wave #it{v}_{2}", 1, pt_min, pt_max + 0.5, 1, 0., 0.8);
    fr->SetXTitle("p_{T} (GeV/#it{c})");
    fr->SetYTitle("#it{v}_{2}");
    fr->Draw();
    gv2.SetMarkerStyle(kFullCircle);
    gv2.SetMarkerSize(0.5);
    gv2.DrawClone("l");

    TCanvas *c2 = new TCanvas("c2");
    c2->SetLogy();
    double xmin, ymin, xmax, ymax;
    giy.ComputeRange(xmin, ymin, xmax, ymax);
    TH2F *fr2 = new TH2F("frame", "blast wave inv. yield", 1, pt_min, pt_max + 0.5, 1, ymin, ymax);
    fr2->SetXTitle("p_{T} (GeV/#it{c})");
    fr2->SetYTitle("inv. yield");
    fr2->Draw();
    giy.SetMarkerStyle(kFullCircle);
    giy.SetMarkerSize(0.5);
    giy.DrawClone("l");

    //
    // test 3: compare to Retiere-Lisa paper results (Fig. 10)
    //
    const double T_RL = 0.1;
    const double Rho0_RL = 0.9;
    const double Rho2_RL = 0.05;
    const double RxOverRy_RL = 11. / 13.;

    const double m_pion = 0.14;
    const double m_proton = 0.938;

    double inv_yield_RL = 0;
    double v2_RL = 0;

    cout << endl << "Test 3: Comparison to Retiere-Lisa paper" << endl;

    // pt = 0.2 GeV, pions
    bw.calc_blastwave_yield_and_v2_fos1(0.2, m_pion, T_RL, Rho0_RL, Rho2_RL, RxOverRy_RL, inv_yield_RL, v2_RL);
    cout << "pt = 0.2 GeV/c, pions, v2 = " << v2_RL << " (value in the Retiere-Lisa paper: 0.018) " << endl;

    // pt = 0.2 GeV, protons
    bw.calc_blastwave_yield_and_v2_fos1(0.2, m_proton, T_RL, Rho0_RL, Rho2_RL, RxOverRy_RL, inv_yield_RL, v2_RL);
    cout << "pt = 0.2 GeV/c, protons, v2 = " << v2_RL << " (value in the Retiere-Lisa paper: 0.004) " << endl;

    // pt = 0.5 GeV, pions
    bw.calc_blastwave_yield_and_v2_fos1(0.5, m_pion, T_RL, Rho0_RL, Rho2_RL, RxOverRy_RL, inv_yield_RL, v2_RL);
    cout << "pt = 0.5 GeV/c, pions, v2 = " << v2_RL << " (value in the Retiere-Lisa paper: 0.062) " << endl;

    // pt = 0.5 GeV, protons
    bw.calc_blastwave_yield_and_v2_fos1(0.5, m_proton, T_RL, Rho0_RL, Rho2_RL, RxOverRy_RL, inv_yield_RL, v2_RL);
    cout << "pt = 0.5 GeV/c, protons, v2 = " << v2_RL << " (value in the Retiere-Lisa paper: 0.025) " << endl;

    // pt = 0.9 GeV, pions
    bw.calc_blastwave_yield_and_v2_fos1(0.9, m_pion, T_RL, Rho0_RL, Rho2_RL, RxOverRy_RL, inv_yield_RL, v2_RL);
    cout << "pt = 0.9 GeV/c, pions, v2 = " << v2_RL << " (value in the Retiere-Lisa paper: 0.112) " << endl;

    // pt = 0.9 GeV, protons
    bw.calc_blastwave_yield_and_v2_fos1(0.9, m_proton, T_RL, Rho0_RL, Rho2_RL, RxOverRy_RL, inv_yield_RL, v2_RL);
    cout << "pt = 0.9 GeV/c, protons, v2 = " << v2_RL << " (value in the Retiere-Lisa paper: 0.067) " << endl;

    //
    // tests for freeze-out hyper-surface 2
    //

    // test 1 boost
    bw.calc_blastwave_yield_and_v2_boost(pt, m, T, rho0, rho2, RxOverRy, inv_yield, v2);

    cout << endl;
    cout << "Test 1 (boost): inv. yield and v2 at a fixed point (for comparison with Mathematica result)" << endl;
    cout << "v2 = " << v2 << endl;
    cout << "invariant yield = " << inv_yield << endl;

    TGraph giy_boost(i_max);
    TGraph gv2_boost(i_max);

    cout << endl << "Test 2 (boost): Performance test (plot v2 vs pt)" << endl;

    TStopwatch t_boost;
    t_boost.Start();

    for (int i = 0; i < i_max; ++i) {
        double pt_BW = pt_min + i * dpt + 0.2 * dpt;
        double inv_yield_BW = 0;
        double v2_BW = 0;
        bw.calc_blastwave_yield_and_v2_boost(pt_BW, m_BW, T_BW, Rho0_BW, Rho2_BW, RxOverRy_BW, inv_yield_BW, v2_BW);
        giy_boost.SetPoint(i, pt_BW, inv_yield_BW);
        gv2_boost.SetPoint(i, pt_BW, v2_BW);
    }

    t_boost.Stop();
    t_boost.Print();

    // plot v2 vs pt
    gStyle->SetOptStat(0);

    TCanvas *c1_boost = new TCanvas("c1_boost");
    gv2_boost.ComputeRange(xmin, ymin, xmax, ymax);
    TH2F *fr_boost = new TH2F("frame", "blast wave #it{v}_{2} (boost)", 1, pt_min, pt_max + 0.5, 1, ymin, ymax);
    fr_boost->SetXTitle("p_{T} (GeV/#it{c})");
    fr_boost->SetYTitle("#it{v}_{2}");
    fr_boost->Draw();
    gv2_boost.SetMarkerStyle(kFullCircle);
    gv2_boost.SetMarkerSize(0.5);
    gv2_boost.DrawClone("l");

    TCanvas *c2_boost = new TCanvas("c2_boost");
    c2_boost->SetLogy();
    // double xmin, ymin, xmax, ymax;
    giy_boost.ComputeRange(xmin, ymin, xmax, ymax);
    TH2F *fr2_boost = new TH2F("frame", "blast wave inv. yield (boost)", 1, pt_min, pt_max + 0.5, 1, ymin, ymax);
    fr2_boost->SetXTitle("p_{T} (GeV/#it{c})");
    fr2_boost->SetYTitle("inv. yield");
    fr2_boost->Draw();
    giy_boost.SetMarkerStyle(kFullCircle);
    giy_boost.SetMarkerSize(0.5);
    giy_boost.DrawClone("l");
}