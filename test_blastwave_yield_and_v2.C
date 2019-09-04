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

    double PhiB = TMath::ATan(RxOverRy * TMath::Tan(PhiHat));  // boost angle
    double rho = rho0 + rho2 * TMath::Cos(2 * PhiB);  // transverse rapidity

	const double A = 1e12;

    return A * rHat * TMath::BesselI(2, (pt * TMath::SinH(rHat * rho)) / T) *
           TMath::BesselK(1, (TMath::Sqrt(TMath::Power(m, 2) + TMath::Power(pt, 2)) * TMath::CosH(rHat * rho)) / T) *
           TMath::Cos(2 * PhiB);
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

	double PhiB = TMath::ATan(RxOverRy * TMath::Tan(PhiHat));  // boost angle
    double rho = rho0 + rho2 * TMath::Cos(2 * PhiB);  // transverse rapidity

    const double A = 1e12; // arbitrary factor to improve numerical stability

    return A * rHat * TMath::BesselI(0, (pt * TMath::SinH(rHat * rho)) / T) *
           TMath::BesselK(1, (TMath::Sqrt(TMath::Power(m, 2) + TMath::Power(pt, 2)) *
                              TMath::CosH(rHat * (rho0 + rho2 * TMath::Cos(2 * PhiB)))) /
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
    ig_num.SetRelTolerance(1e-5);
    
    // wrapper function for numerical integration of v2 denominator
    ROOT::Math::WrappedParamFunction<> w_v2_den(&v2_denominator, 2, 6);
    w_v2_den.SetParameters(pars);
    ROOT::Math::AdaptiveIntegratorMultiDim ig_den;
    ig_den.SetFunction(w_v2_den);
    ig_den.SetRelTolerance(1e-5);
    
    // integration range
    double xmin[2] = {0., 0.};
    double xmax[2] = {1., 2. * TMath::Pi()};

    // integrate
    double v2_num = ig_num.Integral(xmin, xmax);
    // if (ig_num.Status() != 0) cout << ig_num.Status() << endl;
	
    double v2_den = ig_den.Integral(xmin, xmax);
	// if (ig_den.Status() != 0) cout << ig_den.Status() << endl;
	
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

   	blastwave_yield_and_v2(pt, m, T, rho0, rho2, RxOverRy, inv_yield, v2);
 
    cout << "v2 = " << v2 << endl;
    cout << "invariant yield = " << inv_yield << endl;


    //
    // test 2: plot v2 vs pt
    //
    const double T_BW = 0.401757;
    const double Rho0_BW = 1.;
    const double Rho2_BW = 0.167473;
    const double RxOverRy_BW = 0.58366;
    // const double m_BW = 9.460;
    const double m_BW = 3.096;

    const double pt_min = 0;
    const double pt_max = 16.;
    const double dpt = 0.05;

    const int i_max = (pt_max - pt_min) / dpt + 1;
    TGraph gv2(i_max);

	TStopwatch t;
	t.Start(); 

    for (int i=0; i<i_max; ++i) {
    	double pt_BW = pt_min + i * dpt + 0.2 * dpt;
    	double inv_yield_BW = 0;
    	double v2_BW = 0;
    	blastwave_yield_and_v2(pt_BW, m_BW, T_BW, Rho0_BW, Rho2_BW, RxOverRy_BW, inv_yield_BW, v2_BW);
    	gv2.SetPoint(i, pt_BW, v2_BW);
    }

   	t.Stop();
	t.Print(); 

	// plot v2 vs pt
	gStyle->SetOptStat(0);
	TCanvas * c1 = new TCanvas("c1");
	TH2F* fr = new TH2F("frame", "blast wave #it{v}_{2}", 1, pt_min, pt_max + 0.5, 1, 0., 0.8);
	fr->SetXTitle("p_{T} (GeV/#it{c})");
	fr->SetYTitle("#it{v}_{2}");
	fr->Draw();
	gv2.SetMarkerStyle(kCircle);
	gv2.DrawClone("p");

}