#include "functions_BW.h"
#include "Fit/Fitter.h"
#include "Fit/Chi2FCN.h"
#include <Math/Functor.h>
#include "Math/WrappedParamFunction.h"
#include "Math/AdaptiveIntegratorMultiDim.h"

double v2_numerator(const double *x, const double *p);
double v2_denominator(const double *x, const double *p);
void blastwave_yield_and_v2(const double pt, const double m, const double T, const double rho0, const double rho2, const double RxOverRy, double &inv_yield, double &v2);

void BlastWaveMinuitFit()
{
  init_data();

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.15);
  c_play->SetTicks(1,1);
  c_play->SetGrid(0,0);
  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-0.1,9.9);
  for(int i_bin = 0; i_bin < 200; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-100.0);
    h_frame->SetBinError(i_bin+1,1.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.05);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(-0.05,0.5);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("v_{2}");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->Draw("pE");
  PlotLine(0.0,9.9,0.0,0.0,1,2,2);


  Int_t plot_centrality = 4;
  // 30-40%
  vec_graphs[plot_centrality] ->SetMarkerColor(arr_color_mass[0]);
  vec_graphs[plot_centrality] ->Draw("same P"); // pions
  //vec_graphs[plot_centrality+7] ->Draw("same P"); // charged kaons
  vec_graphs[plot_centrality+14] ->SetMarkerColor(arr_color_mass[1]);
  vec_graphs[plot_centrality+14] ->Draw("same P"); // K0s
  vec_graphs[plot_centrality+28] ->SetMarkerColor(arr_color_mass[2]);
  vec_graphs[plot_centrality+28] ->Draw("same P"); // protons


  auto chi2Function = [&](const Double_t *par) 
  {
    //minimisation function computing the sum of squares of residuals
    double chi2 = 0;
    double mass[3] = {0.140,0.498,0.938};
    int plotId[3] = {4,18,32};
    for(int i_pid = 0; i_pid < 3; ++i_pid)
    {
      int fit_id = plotId[i_pid];
      for(int i_point = 0; i_point < vec_graphs[plot_centrality]->GetN(); ++i_point)
      {
	double v2_data = 0.0;
	double pt_data = 0.0;

	vec_graphs[fit_id]->GetPoint(i_point,pt_data,v2_data);
	double v2_err = vec_graphs[plot_centrality]->GetErrorYhigh(i_point);
	// cout << "i_point = " << i_point << ", pt_data = " << pt_data << "v2 = " << v2_data << " +/- " << v2_err << endl;

	if(pt_data < 2.5)
	{
	  double v2_BW = 0;
	  double inv_yield_BW = 0;

	  // blast wave parameters
	  const double pt_BW = pt_data;         // in GeV
	  const double m = mass[i_pid];       // in GeV
	  const double T = par[0];       // fit parameter: Temp in GeV
	  const double rho0 = par[1];   // fit parameter: transverse rapidity
	  const double rho2 = par[2];    // fit parameter: azimuthal modulation of transverse rapidity
	  const double RxOverRy = par[3]; // fit parameter: ratio of the radii Rx and Ry of the freeze-out ellipse in the transverse plane

	  blastwave_yield_and_v2(pt_BW, m, T, rho0, rho2, RxOverRy, inv_yield_BW, v2_BW);
	  // cout << "i_point = " << i_point << ", pt_BW = " << pt_BW << ", v2_BW = " << v2_BW << endl;
	  double diff   = (v2_data - v2_BW)/v2_err;
	  chi2 += diff*diff;
	}
      }
    }
    return chi2;
  };

  // wrap chi2 funciton in a function object for the fit
  // 3 is the number of fit parameters (size of array par)
  ROOT::Math::Functor fcn(chi2Function,4);
  ROOT::Fit::Fitter fitter;
  // ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);

  double pStart[4] = {0.14,0.925,0.05,0.9};
  fitter.SetFCN(fcn, pStart);
   
  fitter.Config().ParSettings(0).SetName("T"); // set parameters name
  fitter.Config().ParSettings(1).SetName("rho0");
  fitter.Config().ParSettings(2).SetName("rho2");
  fitter.Config().ParSettings(3).SetName("RxOverRy");

  // fitter.Config().ParSettings(0).SetLimits(0.01,0.30); // set parameters range
  fitter.Config().ParSettings(0).Fix();
  fitter.Config().ParSettings(1).SetLimits(0.01,20.0);
  fitter.Config().ParSettings(2).SetLimits(0.01,1.0);
  fitter.Config().ParSettings(3).SetLimits(0.01,10.0);

  fitter.Config().ParSettings(0).SetStepSize(0.01); // set parameters minimizer step size
  fitter.Config().ParSettings(1).SetStepSize(0.01);
  fitter.Config().ParSettings(2).SetStepSize(0.01);
  fitter.Config().ParSettings(3).SetStepSize(0.01);

  fitter.Config().SetMinimizer("Minuit","Migrad"); // set minimizer
  // fitter.Config().SetMinimizer("Minuit2","Migrad"); // set minimizer

  // do the fit 
  bool ok = fitter.FitFCN();
  if (!ok) {
    Error("BWFit","BlastWave Fit failed");
  }   

  const ROOT::Fit::FitResult & result = fitter.Result();
  result.Print(std::cout);
  const double *fitpar = result.GetParams();
  double T_BW        = fitpar[0];
  double Rho0_BW     = fitpar[1];
  double Rho2_BW     = fitpar[2];
  double RxOverRy_BW = fitpar[3];
  cout << "T_BW = " << T_BW << ", Rho0_BW = " << Rho0_BW << ", Rho2_BW = " << Rho2_BW << ", RxOverRy_BW = " << RxOverRy_BW << endl;

  /*
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetMaxIterations(10000);  // for GSL
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
  // create funciton wrapper for minmizer
  // a IMultiGenFunction type
  ROOT::Math::Functor fcn(chi2Function,4);
  // starting point
  double pStart[4] = {0.14,0.925,0.05,0.9};
  double step[4] = {0.01,0.01,0.01,0.01};
  min->SetFunction(fcn);
  // Set the free variables to be minimized!
  // min->SetFixedVariable(0,"T",pStart[0]);
  min->SetLimitedVariable(0,"T",pStart[0],step[0],0.01,0.30);
  min->SetLimitedVariable(1,"rho0",pStart[1], step[1],0.01,20.0);
  min->SetLimitedVariable(2,"rho2",pStart[2], step[2],0.01,1.0);
  min->SetLimitedVariable(3,"RxOverRy_BW",pStart[3], step[3],0.01,10.0);
  // do the minimization
  min->Minimize();
  */
}

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

    return rHat * TMath::BesselI(2, (pt * TMath::SinH(rHat * rho)) / T) *
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

    return rHat * TMath::BesselI(0, (pt * TMath::SinH(rHat * rho)) / T) *
           TMath::BesselK(1, (TMath::Sqrt(TMath::Power(m, 2) + TMath::Power(pt, 2)) *
                              TMath::CosH(rHat * (rho0 + rho2 * TMath::Cos(2 * PhiB)))) /
                                 T);
}

void blastwave_yield_and_v2(const double pt, const double m, const double T, const double rho0, const double rho2,
                            const double RxOverRy, double &inv_yield, double &v2) {

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