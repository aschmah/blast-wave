
#include "TString.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGaxis.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TExec.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TPolyLine.h"
#include "TVector3.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TVirtualFitter.h"
#include "Math/MinimizerOptions.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLCamera.h"
#include "TGLPerspectiveCamera.h"
#include "TGFrame.h"
#include "TGLUtil.h"
#include "TGLLightSet.h"
#include "TGLCameraOverlay.h"
#include "TLorentzVector.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include "TF1.h"
#include "TProfile.h"
#include "TArrow.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "Math/WrappedParamFunction.h"
#include "Math/AdaptiveIntegratorMultiDim.h"

static TFile *inputfile_id;
static TFile *inputfile_JPsi;
static TFile *inputfile_D;
static TFile *inputfile_spectra_id;
static TFile *inputfile_spectra_phi;
static TFile *inputfile_spectra_Omega;
static TFile *inputfle_D_dNdpT;
static vector<TString> arr_labels;
static vector<Int_t>   arr_pid;
static vector<TGraphAsymmErrors*> vec_graphs;
static vector< vector<TGraphAsymmErrors*> > vec_tgae_pT_spectra;

 // pi, K+/-, K0s, <K>, p, phi, Lambda, Xi, Omega
static const Int_t N_v2_vs_pt_BW = 9;
static Double_t mt_m0_cut;
static Double_t pt_cut_BW_global = 3.0;
static Int_t flag_v2_BW_use[N_v2_vs_pt_BW]  = {0};
static Int_t flag_v2_BW_plot[N_v2_vs_pt_BW] = {0};
static TGraphAsymmErrors* tgae_v2_stat_BW[N_v2_vs_pt_BW];
static Double_t Mass[N_v2_vs_pt_BW+3] = {0.13957,0.493677,0.497648,0.49566250,0.938272,1.019460,1.115683,1.32131,1.67245,3.096916,1.86962,9.46030};
//static Double_t Mass[N_v2_vs_pt_BW+2] = {0.13957,0.493677,0.497648,0.49566250,0.938272,1.019460,1.115683,1.32131,1.67245,2.5,1.0};
static Double_t arr_pt_low_cut[N_v2_vs_pt_BW+3];
static Double_t arr_pt_high_cut[N_v2_vs_pt_BW+3];
static Int_t arr_color[N_v2_vs_pt_BW+3] = {kBlack,kGreen+1,kRed,kMagenta+1,kCyan+1,kOrange,kRed,kGray,kYellow+2,kRed,kMagenta,kGreen+1};

static const Int_t    N_masses         = 8;
static TH2D* h2D_geometric_shape = NULL;
static TF1 *f_LevyFitFunc        = NULL;
static TF1 *f_FitBessel          = NULL;
static TF1 *f_JetPtFunc          = NULL;
static Double_t arr_quark_mass_meson[N_masses]         = {0.13957,0.497648,0.938272,1.019460,1.67245,1.86962,3.096916,9.46030};
static Int_t    arr_color_mass[N_masses]               = {kBlack,kGreen+1,kRed,kMagenta+1,kCyan+1,kOrange,kYellow+2,kAzure-2};
static const Double_t R_Pb = 5.4946; // fm
static TH2D* h2D_density_Glauber;

static vector<TGraph*> vec_tg_v2_vs_pT_Mathematica;
static vector<TGraph*> vec_tg_dNdpT_vs_pT_Mathematica;
static TGraphErrors* tg_Upsilon_v2_vs_pT;
static TGraphAsymmErrors* tg_JPsi_v2_vs_pT;
static TGraphAsymmErrors* tg_D0_v2_vs_pT;
static vector<TH1F*> h_dN_dpT_mesons;
static vector< vector<TGraphErrors*> > tge_JPsi_spectra;
static TGraphErrors* tge_JPsi_forward_spectrum_stat;
static TGraphErrors* tge_JPsi_forward_spectrum_syst;
static TGraphAsymmErrors* tge_D_dNdpT;
static TGraphAsymmErrors* tge_phi_dNdpT;
static TGraphAsymmErrors* tge_Omega_dNdpT[3];
static TH1D* h_dNdpT_best = NULL;
static vector<TGraphErrors*> vec_tge_v2_vs_pT_560_pid;
static TString label_pid_spectra[8] = {"#pi","K","p","#phi","#Omega","D^{0}","J/#Psi","#Upsilon"};
static TString label_full_pid_spectra[8] = {"pi","K","p","phi","Omega","D0","J/Psi","Upsilon"};
static TString label_v2_dNdpT[2] = {"v2","dNdpT"};

static Double_t Temp_loop_start  = 0.08;
static Double_t rho_0_loop_start = 0.3;
static Double_t rho_a_loop_start = 0.0;
static Double_t Delta_Temp       = 0.02;
static Double_t Delta_rho_0      = 0.125;
static Double_t Delta_rho_a      = 0.05;
static Double_t arr_R_x[9]       = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
static Double_t arr_f_boost[9]   = {0.0,0.05,0.1,0.15,0.2,0.4,0.6,0.8,1.0};
static Int_t    arr_color_line_mass[8] = {kBlack,kGreen,kRed,kMagenta,kCyan,kOrange,kAzure,kGray};

static TGraphAsymmErrors* tgae_v2_vs_pT_mesons_data[8]; // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
static TGraphAsymmErrors* tgae_v2_vs_pT_mesons_data_copy[8]; // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
static TGraphAsymmErrors* tgae_v2_vs_pT_mesons_data_copyB[8]; // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
static TGraphAsymmErrors* tgae_dN_dpT_mesons_data[8];    // pi, K, p, phi, Omega, D0, J/Psi, Upsilon

static Int_t N_calls_BW_ana = 0;

static TMatrixD TL_Matrix_Lambda_L(4,4);
static TMatrixD TL_Matrix_Lambda_LInv(4,4);
static TMatrixD TL_Matrix_Lambda_R(4,4);
static TMatrixD TL_Matrix_Lambda_RInv(4,4);
static TMatrixD TL_Matrix_Lambda_T(4,4);
static TMatrixD TL_Matrix_Lambda_TInv(4,4);


//------------------------------------------------------------------------------------------------------------
static const Float_t Pi = TMath::Pi();
static TRandom ran;
static TString HistName;
static char NoP[50];
//------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------
void Init_Lorentz_Matrices(Double_t y_z, Double_t phis, Double_t y_tp)
{
    //--------------------------------------------
    TL_Matrix_Lambda_L[0][0] = TMath::CosH(y_z);
    TL_Matrix_Lambda_L[0][1] = 0.0;
    TL_Matrix_Lambda_L[0][2] = 0.0;
    TL_Matrix_Lambda_L[0][3] = -TMath::SinH(y_z);

    TL_Matrix_Lambda_L[1][0] = 0.0;
    TL_Matrix_Lambda_L[1][1] = 1.0;
    TL_Matrix_Lambda_L[1][2] = 0.0;
    TL_Matrix_Lambda_L[1][3] = 0.0;

    TL_Matrix_Lambda_L[2][0] = 0.0;
    TL_Matrix_Lambda_L[2][1] = 0.0;
    TL_Matrix_Lambda_L[2][2] = 1.0;
    TL_Matrix_Lambda_L[2][3] = 0.0;

    TL_Matrix_Lambda_L[3][0] = -TMath::SinH(y_z);
    TL_Matrix_Lambda_L[3][1] = 0.0;
    TL_Matrix_Lambda_L[3][2] = 0.0;
    TL_Matrix_Lambda_L[3][3] = TMath::CosH(y_z);
    //--------------------------------------------


    //--------------------------------------------
    TL_Matrix_Lambda_LInv[0][0] = TMath::CosH(-y_z);
    TL_Matrix_Lambda_LInv[0][1] = 0.0;
    TL_Matrix_Lambda_LInv[0][2] = 0.0;
    TL_Matrix_Lambda_LInv[0][3] = -TMath::SinH(-y_z);

    TL_Matrix_Lambda_LInv[1][0] = 0.0;
    TL_Matrix_Lambda_LInv[1][1] = 1.0;
    TL_Matrix_Lambda_LInv[1][2] = 0.0;
    TL_Matrix_Lambda_LInv[1][3] = 0.0;

    TL_Matrix_Lambda_LInv[2][0] = 0.0;
    TL_Matrix_Lambda_LInv[2][1] = 0.0;
    TL_Matrix_Lambda_LInv[2][2] = 1.0;
    TL_Matrix_Lambda_LInv[2][3] = 0.0;

    TL_Matrix_Lambda_LInv[3][0] = -TMath::SinH(-y_z);
    TL_Matrix_Lambda_LInv[3][1] = 0.0;
    TL_Matrix_Lambda_LInv[3][2] = 0.0;
    TL_Matrix_Lambda_LInv[3][3] = TMath::CosH(-y_z);
    //--------------------------------------------



    //--------------------------------------------
    TL_Matrix_Lambda_R[0][0] = 1.0;
    TL_Matrix_Lambda_R[0][1] = 0.0;
    TL_Matrix_Lambda_R[0][2] = 0.0;
    TL_Matrix_Lambda_R[0][3] = 0.0;

    TL_Matrix_Lambda_R[1][0] = 0.0;
    TL_Matrix_Lambda_R[1][1] = TMath::Cos(phis);
    TL_Matrix_Lambda_R[1][2] = TMath::Sin(phis);
    TL_Matrix_Lambda_R[1][3] = 0.0;

    TL_Matrix_Lambda_R[2][0] = 0.0;
    TL_Matrix_Lambda_R[2][1] = -TMath::Sin(phis);
    TL_Matrix_Lambda_R[2][2] = TMath::Cos(phis);
    TL_Matrix_Lambda_R[2][3] = 0.0;

    TL_Matrix_Lambda_R[3][0] = 0.0;
    TL_Matrix_Lambda_R[3][1] = 0.0;
    TL_Matrix_Lambda_R[3][2] = 0.0;
    TL_Matrix_Lambda_R[3][3] = 1.0;
    //--------------------------------------------



    //--------------------------------------------
    TL_Matrix_Lambda_RInv[0][0] = 1.0;
    TL_Matrix_Lambda_RInv[0][1] = 0.0;
    TL_Matrix_Lambda_RInv[0][2] = 0.0;
    TL_Matrix_Lambda_RInv[0][3] = 0.0;

    TL_Matrix_Lambda_RInv[1][0] = 0.0;
    TL_Matrix_Lambda_RInv[1][1] = TMath::Cos(-phis);
    TL_Matrix_Lambda_RInv[1][2] = TMath::Sin(-phis);
    TL_Matrix_Lambda_RInv[1][3] = 0.0;

    TL_Matrix_Lambda_RInv[2][0] = 0.0;
    TL_Matrix_Lambda_RInv[2][1] = -TMath::Sin(-phis);
    TL_Matrix_Lambda_RInv[2][2] = TMath::Cos(-phis);
    TL_Matrix_Lambda_RInv[2][3] = 0.0;

    TL_Matrix_Lambda_RInv[3][0] = 0.0;
    TL_Matrix_Lambda_RInv[3][1] = 0.0;
    TL_Matrix_Lambda_RInv[3][2] = 0.0;
    TL_Matrix_Lambda_RInv[3][3] = 1.0;
    //--------------------------------------------



    //--------------------------------------------
    TL_Matrix_Lambda_T[0][0] = TMath::CosH(y_tp);
    TL_Matrix_Lambda_T[0][1] = -TMath::SinH(y_tp);
    TL_Matrix_Lambda_T[0][2] = 0.0;
    TL_Matrix_Lambda_T[0][3] = 0.0;

    TL_Matrix_Lambda_T[1][0] = -TMath::SinH(y_tp);
    TL_Matrix_Lambda_T[1][1] = TMath::CosH(y_tp);
    TL_Matrix_Lambda_T[1][2] = 0.0;
    TL_Matrix_Lambda_T[1][3] = 0.0;

    TL_Matrix_Lambda_T[2][0] = 0.0;
    TL_Matrix_Lambda_T[2][1] = 0.0;
    TL_Matrix_Lambda_T[2][2] = 1.0;
    TL_Matrix_Lambda_T[2][3] = 0.0;

    TL_Matrix_Lambda_T[3][0] = 0.0;
    TL_Matrix_Lambda_T[3][1] = 0.0;
    TL_Matrix_Lambda_T[3][2] = 0.0;
    TL_Matrix_Lambda_T[3][3] = 1.0;
    //--------------------------------------------


    //--------------------------------------------
    TL_Matrix_Lambda_TInv[0][0] = TMath::CosH(-y_tp);
    TL_Matrix_Lambda_TInv[0][1] = -TMath::SinH(-y_tp);
    TL_Matrix_Lambda_TInv[0][2] = 0.0;
    TL_Matrix_Lambda_TInv[0][3] = 0.0;

    TL_Matrix_Lambda_TInv[1][0] = -TMath::SinH(-y_tp);
    TL_Matrix_Lambda_TInv[1][1] = TMath::CosH(-y_tp);
    TL_Matrix_Lambda_TInv[1][2] = 0.0;
    TL_Matrix_Lambda_TInv[1][3] = 0.0;

    TL_Matrix_Lambda_TInv[2][0] = 0.0;
    TL_Matrix_Lambda_TInv[2][1] = 0.0;
    TL_Matrix_Lambda_TInv[2][2] = 1.0;
    TL_Matrix_Lambda_TInv[2][3] = 0.0;

    TL_Matrix_Lambda_TInv[3][0] = 0.0;
    TL_Matrix_Lambda_TInv[3][1] = 0.0;
    TL_Matrix_Lambda_TInv[3][2] = 0.0;
    TL_Matrix_Lambda_TInv[3][3] = 1.0;
    //--------------------------------------------
}
//------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------
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
    double A = p[6]; // arbitrary factor, can be adjusted to improve numerical stability

    double PhiB = TMath::ATan(RxOverRy * TMath::Tan(PhiHat)); // boost angle
    double rho = rho0 + rho2 * TMath::Cos(2 * PhiB);          // transverse rapidity

    return A * rHat * TMath::BesselI(2, (pt * TMath::SinH(rHat * rho)) / T) *
           TMath::BesselK(1, (TMath::Sqrt(TMath::Power(m, 2) + TMath::Power(pt, 2)) * TMath::CosH(rHat * rho)) / T) *
           TMath::Cos(2 * PhiB);
}
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
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
    double A = p[6]; // arbitrary factor, can be adjusted to improve numerical stability

    double PhiB = TMath::ATan(RxOverRy * TMath::Tan(PhiHat)); // boost angle
    double rho = rho0 + rho2 * TMath::Cos(2 * PhiB);          // transverse rapidity

    return A * rHat * TMath::BesselI(0, (pt * TMath::SinH(rHat * rho)) / T) *
           TMath::BesselK(1, (TMath::Sqrt(TMath::Power(m, 2) + TMath::Power(pt, 2)) *
                              TMath::CosH(rHat * (rho0 + rho2 * TMath::Cos(2 * PhiB)))) /
                                 T);
}
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
void blastwave_yield_and_v2(const double &pt, const double &m, const double &T, const double &rho0, const double &rho2,
                            const double &RxOverRy, double &inv_yield, double &v2) {

    // mass and T dependent scale factor prop. to total yield used to improve numerical stability
    // of the invariant yield
    double sf = T * m * m * TMath::BesselK(2, m / T);

    // blast wave parameters:
    // the last number (par[6]) is an arbitrary normalization which we will adjust
    // in order to have a good numerical stability for different masses and pt
    double pars[7] = {pt, m, T, rho0, rho2, RxOverRy, 1. / sf};

    // wrapper functions for v2 numerator and denominator
    ROOT::Math::WrappedParamFunction<> w_v2_num(&v2_numerator, 2, 7);
    ROOT::Math::WrappedParamFunction<> w_v2_den(&v2_denominator, 2, 7);

    // define integrator
    ROOT::Math::AdaptiveIntegratorMultiDim ig;
    ig.SetRelTolerance(1e-6);

    // integration range
    double xmin[2] = {0., 0.};
    double xmax[2] = {1., 2. * TMath::Pi()};

    // We calculate the invariant yield 1/(2 pi pt) dN/(dpt dy) (with arbitrary constant pre-factor)
    // and determine the optimal normalization for numerical stability in the v2 calculation
    w_v2_den.SetParameters(pars);
    ig.SetFunction(w_v2_den);
    inv_yield = sf * TMath::Sqrt(m * m + pt * pt) * ig.Integral(xmin, xmax);
    pars[6] = 1. / inv_yield;
    w_v2_num.SetParameters(pars);
    w_v2_den.SetParameters(pars);

    // integrate
    ig.SetFunction(w_v2_num);
    double v2_num = ig.Integral(xmin, xmax);
    // if (ig_num.Status() != 0) cout << ig_num.Status() << endl;

    ig.SetFunction(w_v2_den);
    double v2_den = ig.Integral(xmin, xmax);
    // if (ig_den.Status() != 0) cout << ig_den.Status() << endl;

    // cout << pt << " " << v2_den << endl;
    // cout << pt << " " << inv_yield << endl;

    if (v2_den != 0) {
        v2 = v2_num / v2_den;
    } else {
        cout << "WARNING: v2 denominator zero!!!" << endl;
    }
}
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
void Multiply_pT_norm_to_integral_tgae(TGraphAsymmErrors* tgae_in = NULL, Int_t flag_mult_pT = 1)
{
    // Modify a TGraphAsymmErrors by multiplying each value with pT
    // In a second step the spectrum is normalized to its integral
    // Multiply with pt
    Double_t integral = 0.0;
    for(Int_t i_point = 0; i_point < tgae_in ->GetN(); i_point++)
    {
        Double_t pT,y_val,err_X_high,err_X_low,err_Y_high,err_Y_low;
        tgae_in ->GetPoint(i_point,pT,y_val);
        err_X_high = tgae_in ->GetErrorXhigh(i_point);
        err_X_low  = tgae_in ->GetErrorXlow(i_point);
        err_Y_high = tgae_in ->GetErrorYhigh(i_point);
        err_Y_low  = tgae_in ->GetErrorYlow(i_point);

        if(flag_mult_pT)
        {
            tgae_in ->SetPoint(i_point,pT,y_val*pT);
            tgae_in ->SetPointEYhigh(i_point,err_Y_high*pT);
            tgae_in ->SetPointEYlow(i_point,err_Y_low*pT);
            integral += y_val*pT*(err_X_low + err_X_high); // bin content * bin width
        }
        else
        {
            tgae_in ->SetPoint(i_point,pT,y_val);
            tgae_in ->SetPointEYhigh(i_point,err_Y_high);
            tgae_in ->SetPointEYlow(i_point,err_Y_low);
            integral += y_val*(err_X_low + err_X_high); // bin content * bin width
        }

    }

    if(integral >= 0.0)
    {
        // Normalize to integral
        for(Int_t i_point = 0; i_point < tgae_in ->GetN(); i_point++)
        {
            Double_t pT,y_val,err_X_high,err_X_low,err_Y_high,err_Y_low;
            tgae_in ->GetPoint(i_point,pT,y_val);
            err_X_high = tgae_in ->GetErrorXhigh(i_point);
            err_X_low  = tgae_in ->GetErrorXlow(i_point);
            err_Y_high = tgae_in ->GetErrorYhigh(i_point);
            err_Y_low  = tgae_in ->GetErrorYlow(i_point);

            tgae_in ->SetPoint(i_point,pT,y_val/integral);
            tgae_in ->SetPointEYhigh(i_point,err_Y_high/integral);
            tgae_in ->SetPointEYlow(i_point,err_Y_low/integral);
        }
    }
}
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
TGraphAsymmErrors* Add_tgae(TGraphAsymmErrors* tgae_inA = NULL, TGraphAsymmErrors* tgae_inB = NULL)
{
    // Add two TGraphAsymmErrors with identical number of points and positions

    TGraphAsymmErrors* tgae_inAB[2] = {tgae_inA,tgae_inB};
    TGraphAsymmErrors* tgae_out = new TGraphAsymmErrors();

    for(Int_t i_point = 0; i_point < tgae_inAB[0] ->GetN(); i_point++)
    {
        Double_t pT[2],y_val[2],err_X_high[2],err_X_low[2],err_Y_high[2],err_Y_low[2];
        for(Int_t in_AB = 0; in_AB < 2; in_AB++)
        {
            tgae_inAB[in_AB] ->GetPoint(i_point,pT[in_AB],y_val[in_AB]);
            err_X_high[in_AB] = tgae_inAB[in_AB] ->GetErrorXhigh(i_point);
            err_X_low[in_AB]  = tgae_inAB[in_AB] ->GetErrorXlow(i_point);
            err_Y_high[in_AB] = tgae_inAB[in_AB] ->GetErrorYhigh(i_point);
            err_Y_low[in_AB]  = tgae_inAB[in_AB] ->GetErrorYlow(i_point);
        }

        Double_t new_y_val = y_val[0]+y_val[1];
        tgae_out ->SetPoint(i_point,pT[0],new_y_val);
        Double_t err_Y_high_out = TMath::Sqrt(TMath::Power(err_Y_high[0],2.0) + TMath::Power(err_Y_high[1],2.0));
        Double_t err_Y_low_out = TMath::Sqrt(TMath::Power(err_Y_low[0],2.0) + TMath::Power(err_Y_low[1],2.0));
        tgae_out ->SetPointEYhigh(i_point,err_Y_high_out);
        tgae_out ->SetPointEYlow(i_point,err_Y_low_out);
        tgae_out ->SetPointEXhigh(i_point,err_X_high[0]);
        tgae_out ->SetPointEXlow(i_point,err_X_low[0]);

        //printf("i_point: %d, pT: %4.3f, new_y_val: %4.3f, y_vals: {%4.3f, %4.3f} \n",i_point,pT[0],new_y_val,y_val[0],y_val[1]);
    }

    return tgae_out;
}
//------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PtFitBessel(Double_t* x_val, Double_t* par)
{
    Double_t x, y, m0, Temp, Ampl, shift;
    m0    = par[0];
    Temp  = par[1];
    Ampl  = par[2];
    shift = par[3];
    x = x_val[0];
    y = Ampl*x*sqrt(x*x+m0*m0)*TMath::BesselK1(sqrt(x*x+m0*m0)/Temp);
    return y;
}
//----------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
void init_pT_spectra_data()
{
    f_FitBessel = new TF1("f_FitBessel",PtFitBessel,0.0,10.0,4);

    // https://www.hepdata.net/record/ins1377750
    printf("Initialize pT spectra data \n");

    inputfle_D_dNdpT        = TFile::Open("./Data/HEPData-ins1394580-v1-root.root"); // https://www.hepdata.net/record/ins1394580;
    inputfile_spectra_id    = TFile::Open("./Data/HEP_ALICE_PID_pT_spectra.root"); // pi, K, p
    inputfile_spectra_phi   = TFile::Open("./Data/HEPData-ins1511864-v1-root_phi_KStar_dNdpT_2.76TeV.root"); // phi
    inputfile_spectra_Omega = TFile::Open("./Data/HEPData-ins1243865-v1-root_omegas_dNdpT_2.76TeV.root"); // Omega


    tge_D_dNdpT        = (TGraphAsymmErrors*)inputfle_D_dNdpT        ->Get(Form("Table %d/Graph1D_y%d",4,1)); // 30-50%
    tge_phi_dNdpT      = (TGraphAsymmErrors*)inputfile_spectra_phi   ->Get(Form("Table %d/Graph1D_y%d",13,1)); // 30-40%
    tge_Omega_dNdpT[0] = (TGraphAsymmErrors*)inputfile_spectra_Omega ->Get(Form("Table %d/Graph1D_y%d",8,1)); // 20-40%, Omega-
    tge_Omega_dNdpT[1] = (TGraphAsymmErrors*)inputfile_spectra_Omega ->Get(Form("Table %d/Graph1D_y%d",8,2)); // 20-40%, Omega+

    tge_Omega_dNdpT[2] = Add_tgae(tge_Omega_dNdpT[0],tge_Omega_dNdpT[1]);

    Multiply_pT_norm_to_integral_tgae(tge_Omega_dNdpT[2],0);

    const Int_t arr_centrality_low[11]        = {0,5,10,20,40,60,20,30,40,50,5};
    const Int_t arr_centrality_high[11]       = {5,10,20,40,60,80,30,40,50,60,60};
    TString label_pid[3] = {"pi","K","p"};

    vec_tgae_pT_spectra.resize(3); // pi, K, p
    for(Int_t i_pid = 0; i_pid < 3; i_pid++)
    {
        TString label = label_pid[i_pid];
        // 0-5, 5-10, 10-20, 20-40, 40-60, 60-80, 20-30, 30-40, 40-50, (50-60), 5-60
        for(Int_t i = 0; i < 6; i++)
        {
            label += Form("_%d_%d",arr_centrality_low[i],arr_centrality_high[i]);
            vec_tgae_pT_spectra[i_pid].push_back((TGraphAsymmErrors*)inputfile_spectra_id->Get(Form("Table %d/Graph1D_y%d",i_pid+1,i+1)));
            vec_tgae_pT_spectra[i_pid][i]->SetName(label.Data());
            vec_tgae_pT_spectra[i_pid][i]->SetMarkerColor(arr_color[i]);
            vec_tgae_pT_spectra[i_pid][i]->SetMarkerStyle(20);
            vec_tgae_pT_spectra[i_pid][i]->SetMarkerSize(0.8);
            vec_tgae_pT_spectra[i_pid][i]->GetXaxis()->CenterTitle();
            vec_tgae_pT_spectra[i_pid][i]->GetYaxis()->CenterTitle();
            vec_tgae_pT_spectra[i_pid][i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            vec_tgae_pT_spectra[i_pid][i]->GetYaxis()->SetTitle("1/p_{T} dN/dp_{T} (GeV/c)^{-2}");
            vec_tgae_pT_spectra[i_pid][i]->SetTitle("");
            vec_tgae_pT_spectra[i_pid][i]->GetXaxis()->SetTitleOffset(1.2);
            vec_tgae_pT_spectra[i_pid][i]->GetYaxis()->SetTitleOffset(1.2);
            vec_tgae_pT_spectra[i_pid][i]->GetXaxis()->SetLabelSize(0.06);
            vec_tgae_pT_spectra[i_pid][i]->GetYaxis()->SetLabelSize(0.06);
            vec_tgae_pT_spectra[i_pid][i]->GetXaxis()->SetTitleSize(0.06);
            vec_tgae_pT_spectra[i_pid][i]->GetYaxis()->SetTitleSize(0.06);
            vec_tgae_pT_spectra[i_pid][i]->GetXaxis()->SetNdivisions(505,'N');
            vec_tgae_pT_spectra[i_pid][i]->GetYaxis()->SetNdivisions(505,'N');
        }
    }

    for(Int_t i_pid = 0; i_pid < 3; i_pid++)
    {
        TString label = label_pid[i_pid];
        // 0-5, 5-10, 10-20, 20-40, 40-60, 60-80, 20-30, 30-40, 40-50, (50-60)
        for(Int_t i = 6; i < 9; i++)
        {
            label += Form("_%d_%d",arr_centrality_low[i],arr_centrality_high[i]);
            vec_tgae_pT_spectra[i_pid].push_back((TGraphAsymmErrors*)inputfile_spectra_id->Get(Form("Table %d/Graph1D_y%d",i_pid+4,i-5)));
            vec_tgae_pT_spectra[i_pid][i]->SetName(label.Data());
            vec_tgae_pT_spectra[i_pid][i]->SetMarkerColor(arr_color[i]);
            vec_tgae_pT_spectra[i_pid][i]->SetMarkerStyle(20);
            vec_tgae_pT_spectra[i_pid][i]->SetMarkerSize(0.8);
            vec_tgae_pT_spectra[i_pid][i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            vec_tgae_pT_spectra[i_pid][i]->GetYaxis()->SetTitle("1/p_{T} dN/dp_{T} (GeV/c)^{-2}");
        }
    }

    // Make 50-60%
    for(Int_t i_pid = 0; i_pid < 3; i_pid++)
    {
        TString label = label_pid[i_pid];
        // 0-5, 5-10, 10-20, 20-40, 40-60, 60-80, 20-30, 30-40, 40-50, (50-60)
        Int_t i = 9; // 50-60% = 40-60% - 40-50%
        Int_t i_40_60 = 4; // 40-60%
        Int_t i_40_50 = 8; // 40-50%
        label += Form("_%d_%d",arr_centrality_low[i],arr_centrality_high[i]);
        vec_tgae_pT_spectra[i_pid].push_back((TGraphAsymmErrors*)vec_tgae_pT_spectra[i_pid][i_40_60]);
        vec_tgae_pT_spectra[i_pid][i]->SetName(label.Data());
        vec_tgae_pT_spectra[i_pid][i]->SetMarkerColor(arr_color[i]);
        vec_tgae_pT_spectra[i_pid][i]->SetMarkerStyle(20);
        vec_tgae_pT_spectra[i_pid][i]->SetMarkerSize(0.8);
        vec_tgae_pT_spectra[i_pid][i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        vec_tgae_pT_spectra[i_pid][i]->GetYaxis()->SetTitle("1/p_{T} dN/dp_{T} (GeV/c)^{-2}");

        for(Int_t i_pT = 0; i_pT < vec_tgae_pT_spectra[i_pid][i]->GetN(); i_pT++)
        {
            Double_t x_val, y_val;
            vec_tgae_pT_spectra[i_pid][i]       ->GetPoint(i_pT,x_val,y_val);
            y_val *= 2.0; // 40-60% was divided by the number of events in a 20% bin, relative to the 10% one in 40-50%
            Double_t y_val_40_50 = vec_tgae_pT_spectra[i_pid][i_40_50] ->Eval(x_val);

            vec_tgae_pT_spectra[i_pid][i]       ->SetPoint(i_pT,x_val,y_val - y_val_40_50);
            //if(i_pid == 0) printf("i_pT: %d, pT: %4.3f, y_val(40-60%): %4.6f, y_val(40-50%): %4.6f \n",i_pT,x_val,y_val,y_val_40_50);
        }
    }

    // Make 5-60%
    for(Int_t i_pid = 0; i_pid < 3; i_pid++)
    {
        TString label = label_pid[i_pid];
        // 0-5, 5-10, 10-20, 20-40, 40-60, 60-80, 20-30, 30-40, 40-50, (50-60), (5-60)
        Int_t i = 10;
        label += Form("_%d_%d",arr_centrality_low[i],arr_centrality_high[i]);
        vec_tgae_pT_spectra[i_pid].push_back((TGraphAsymmErrors*)vec_tgae_pT_spectra[i_pid][1]);
        vec_tgae_pT_spectra[i_pid][i]->SetName(label.Data());
        vec_tgae_pT_spectra[i_pid][i]->SetMarkerColor(arr_color[i]);
        vec_tgae_pT_spectra[i_pid][i]->SetMarkerStyle(20);
        vec_tgae_pT_spectra[i_pid][i]->SetMarkerSize(0.8);
        vec_tgae_pT_spectra[i_pid][i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        vec_tgae_pT_spectra[i_pid][i]->GetYaxis()->SetTitle("1/p_{T} dN/dp_{T} (GeV/c)^{-2}");

        Double_t arr_scale_cent_width[5] = {1.0,1.0,2.0,4.0,4.0}; // width: 5%, 5%, 10%, 20%, 20%

        for(Int_t i_pT = 0; i_pT < vec_tgae_pT_spectra[i_pid][i]->GetN(); i_pT++)
        {
            Double_t x_val, y_val;
            vec_tgae_pT_spectra[i_pid][i]       ->GetPoint(i_pT,x_val,y_val);
            for(Int_t i_cent = 2; i_cent <= 4; i_cent++)
            {
                Double_t y_val_cent = vec_tgae_pT_spectra[i_pid][i_cent] ->Eval(x_val);
                y_val += y_val_cent*arr_scale_cent_width[i_cent];
            }
            vec_tgae_pT_spectra[i_pid][i]       ->SetPoint(i_pT,x_val,y_val);
        }
    }


    //----------------------------------------------------------------------
    for(Int_t i_pid = 0; i_pid < 3; i_pid++)
    {
        for(Int_t i = 0; i < 11; i++)
        {
            Multiply_pT_norm_to_integral_tgae(vec_tgae_pT_spectra[i_pid][i],1);
        }
    }
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    // phi mesons 2.76 TeV
    Multiply_pT_norm_to_integral_tgae(tge_phi_dNdpT,1);
    //----------------------------------------------------------------------



    //----------------------------------------------------------------------
    // J/Psi
    // Private communication Markus Koehler
    // dN/dpt

    // Mid-rapidity
    tge_JPsi_spectra.resize(3); // 0-20%, 20-40%, 40-90%
    tge_JPsi_spectra[0].resize(2); // stat, syst
    tge_JPsi_spectra[1].resize(2); // stat, syst
    tge_JPsi_spectra[2].resize(2); // stat, syst

    if(1)
    {
        const int nPoints = 4;
        double xval[nPoints] = {0.725,2.15,4.0,7.5};
        double xerr[nPoints] = {0.575,0.85,1.0,2.5};

        double yval[nPoints] = {0.0234115,0.0329794,0.00913332,0.000647037};
        double yerr[nPoints] = {0.00460038,0.00354616,0.00107112,0.000120726};

        TGraphErrors *gElePtStat = new TGraphErrors(nPoints, xval, yval, xerr, yerr);
        tge_JPsi_spectra[0][0] = (TGraphErrors*)gElePtStat ->Clone("tge_JPsi_spectra_0_20_stat");
    }

    if(1)
    {
        const int nPoints = 4;
        double xval[nPoints] = {0.725,2.15,4.0,7.5};
        double xerr[nPoints] = {0.575,0.85,1.0,2.5};

        double yval[nPoints] = {0.0234115,0.0329794,0.00913332,0.000647037};
        //  double yerr[nPoints] = {0.0042887,0.00346681,0.000843378,7.34505e-05};

        double yerr[nPoints] = {TMath::Sqrt(0.090*0.090 + 0.063*0.063 + 0.052*0.052 + 0.015*0.015)*yval[0],
        TMath::Sqrt(0.044*0.044 + 0.045*0.045 + 0.052*0.052 + 0.010*0.010)*yval[1],
        TMath::Sqrt(0.043*0.043 + 0.018*0.018 + 0.065*0.065 + 0.015*0.015)*yval[2],
        TMath::Sqrt(0.043*0.043 + 0.013*0.013 + 0.053*0.053 + 0.020*0.020)*yval[3]};// values from note Table6 and Table7 00-90% reduce fluctuations

        TGraphErrors *gElePtSyst = new TGraphErrors(nPoints, xval, yval, xerr, yerr);
        tge_JPsi_spectra[0][1] = (TGraphErrors*)gElePtSyst ->Clone("tge_JPsi_spectra_0_20_syst");
    }

    if(1)
    {
        const int nPoints = 4;
        double xval[nPoints] = {0.725,2.15,4.0,7.5};
        double xerr[nPoints] = {0.575,0.85,1.0,2.5};

        double yval[nPoints] = {0.00705774,0.0102737,0.00279403,0.000404948};//
        double yerr[nPoints] = {0.00172881,0.00131777,0.000424049,5.67846e-05};//

        TGraphErrors *gElePtStat = new TGraphErrors(nPoints, xval, yval, xerr, yerr);
        tge_JPsi_spectra[1][0] = (TGraphErrors*)gElePtStat ->Clone("tge_JPsi_spectra_20_40_stat");
    }

    if(1)
    {
        const int nPoints = 4;
        double xval[nPoints] = {0.725,2.15,4.0,7.5};
        double xerr[nPoints] = {0.575,0.85,1.0,2.5};

        double yval[nPoints] = {0.00705774,0.0102737,0.00279403,0.000404948};//
        //  double yerr[nPoints] = {0.00134022,0.00115116,0.000394142,2.64464e-05};//


        double yerr[nPoints] = {TMath::Sqrt(0.090*0.090 + 0.063*0.063 + 0.052*0.052 + 0.015*0.015)*yval[0],
        TMath::Sqrt(0.044*0.044 + 0.045*0.045 + 0.052*0.052 + 0.010*0.010)*yval[1],
        TMath::Sqrt(0.043*0.043 + 0.018*0.018 + 0.065*0.065 + 0.015*0.015)*yval[2],
        TMath::Sqrt(0.043*0.043 + 0.013*0.013 + 0.053*0.053 + 0.020*0.020)*yval[3]};// values from note Table6 and Table7 00-90% reduce fluctuations

        TGraphErrors *gElePtSyst = new TGraphErrors(nPoints, xval, yval, xerr, yerr);
        tge_JPsi_spectra[1][1] = (TGraphErrors*)gElePtSyst ->Clone("tge_JPsi_spectra_20_40_syst");
    }

    if(1)
    {
        const int nPoints = 4;
        double xval[nPoints] = {0.725,2.15,4.0,7.5};
        double xerr[nPoints] = {0.575,0.85,1.0,2.5};


        double yval[nPoints] = {0.00138902,0.000802937,0.000570093,7.41637e-05};//
        double yerr[nPoints] = {0.000252756,0.000191895,6.82242e-05,1.13561e-05};//

        TGraphErrors *gElePtStat = new TGraphErrors(nPoints, xval, yval, xerr, yerr);
        tge_JPsi_spectra[2][0] = (TGraphErrors*)gElePtStat ->Clone("tge_JPsi_spectra_40_90_stat");
    }

    if(1)
    {
        const int nPoints = 4;
        double xval[nPoints] = {0.725,2.15,4.0,7.5};
        double xerr[nPoints] = {0.575,0.85,1.0,2.5};

        double yval[nPoints] = {0.00138902,0.000802937,0.000570093,7.41637e-05};
        //  double yerr[nPoints] = {0.000233462,0.000203163,5.98407e-05,5.25889e-06};


        double yerr[nPoints] = {TMath::Sqrt(0.090*0.090 + 0.063*0.063 + 0.052*0.052 + 0.015*0.015)*yval[0],
        TMath::Sqrt(0.044*0.044 + 0.045*0.045 + 0.052*0.052 + 0.010*0.010)*yval[1],
        TMath::Sqrt(0.043*0.043 + 0.018*0.018 + 0.065*0.065 + 0.015*0.015)*yval[2],
        TMath::Sqrt(0.043*0.043 + 0.013*0.013 + 0.053*0.053 + 0.020*0.020)*yval[3]};// values from note Table6 and Table7 00-90% reduce fluctuations

        TGraphErrors *gElePtSyst = new TGraphErrors(nPoints, xval, yval, xerr, yerr);
        tge_JPsi_spectra[2][1] = (TGraphErrors*)gElePtSyst ->Clone("tge_JPsi_spectra_40_90_syst");
    }

    for(Int_t iCent = 0; iCent < 3; iCent++)
    {
        Double_t integral = 0.0;
        for(Int_t iPoint = 0; iPoint < tge_JPsi_spectra[iCent][0] ->GetN(); iPoint++)
        {
            Double_t pT, y_val;
            Double_t err_X_low  = tge_JPsi_spectra[iCent][0] ->GetErrorXlow(iPoint);
            Double_t err_X_high = tge_JPsi_spectra[iCent][0] ->GetErrorXhigh(iPoint);
            tge_JPsi_spectra[iCent][0] ->GetPoint(iPoint,pT,y_val);

            integral += y_val*(err_X_low + err_X_high);
        }

        if(integral <= 0.0) continue;
        for(Int_t iPoint = 0; iPoint < tge_JPsi_spectra[iCent][0] ->GetN(); iPoint++)
        {
            Double_t pT, y_val;
            tge_JPsi_spectra[iCent][0] ->GetPoint(iPoint,pT,y_val);
            Double_t err_X_low       = tge_JPsi_spectra[iCent][0] ->GetErrorXlow(iPoint);
            Double_t err_X_high      = tge_JPsi_spectra[iCent][0] ->GetErrorXhigh(iPoint);
            Double_t err_Y_low       = tge_JPsi_spectra[iCent][0] ->GetErrorYlow(iPoint);
            Double_t err_Y_high      = tge_JPsi_spectra[iCent][0] ->GetErrorYhigh(iPoint);
            Double_t err_Y_low_syst  = tge_JPsi_spectra[iCent][1] ->GetErrorYlow(iPoint);
            Double_t err_Y_high_syst = tge_JPsi_spectra[iCent][1] ->GetErrorYhigh(iPoint);
            tge_JPsi_spectra[iCent][0] ->SetPoint(iPoint,pT,y_val/integral);
            tge_JPsi_spectra[iCent][1] ->SetPoint(iPoint,pT,y_val/integral);
            tge_JPsi_spectra[iCent][0] ->SetPointError(iPoint,err_X_high,err_Y_high/integral);
            tge_JPsi_spectra[iCent][1] ->SetPointError(iPoint,err_X_high,err_Y_high_syst/integral);
        }
    }



    // Forward rapidity
    //=================================================
    //    Yields values in pt bins 0-20% (2015) 0<pt<12
    //    =================================================
    //    pt= 0-1 Y =   9.19406  +- 0.32491 (stat) +- 0.59610 (syst) +0.13822 (global)
    //    pt= 1-2 Y =   16.08296 +- 0.36600 (stat) +- 0.99704 (syst) +0.24178 (global)
    //    pt= 2-3 Y =   10.93737 +- 0.26496 (stat) +- 0.59763 (syst) +0.16442 (global)
    //    pt= 3-4 Y =   4.96946  +- 0.12741 (stat) +- 0.24519 (syst) +0.07471 (global)
    //    pt= 4-5 Y =   2.01827  +- 0.06734 (stat) +- 0.09296 (syst) +0.03034 (global)
    //    pt= 5-6 Y =   0.89579  +- 0.03562 (stat) +- 0.03838 (syst) +0.01347 (global)
    //    pt= 6-7 Y =   0.37794  +- 0.01633 (stat) +- 0.01744 (syst) +0.00568 (global)
    //    pt= 7-8 Y =   0.16213  +- 0.01158 (stat) +- 0.00721 (syst) +0.00244 (global)
    //    pt = 8-9 Y =  0.09443  +- 0.00714 (stat) +- 0.00465 (syst) +0.00142 (global)
    //    pt= 9-10 Y =  0.04853  +- 0.00528 (stat) +- 0.00236 (syst) +0.00073 (global)
    //    pt= 10-12 Y = 0.02248  +- 0.00187 (stat) +- 0.00115 (syst) +0.00034 (global)

    tge_JPsi_forward_spectrum_stat = new TGraphErrors();
    tge_JPsi_forward_spectrum_stat ->SetName("tge_JPsi_forward_spectrum_stat");
    tge_JPsi_forward_spectrum_stat ->SetPoint(0,0.5,9.19406);
    tge_JPsi_forward_spectrum_stat ->SetPoint(1,1.5,16.0829);
    tge_JPsi_forward_spectrum_stat ->SetPoint(2,2.5,10.9373);
    tge_JPsi_forward_spectrum_stat ->SetPoint(3,3.5,4.96946);
    tge_JPsi_forward_spectrum_stat ->SetPoint(4,4.5,2.01827);
    tge_JPsi_forward_spectrum_stat ->SetPoint(5,5.5,0.89579);
    tge_JPsi_forward_spectrum_stat ->SetPoint(6,6.5,0.37794);
    tge_JPsi_forward_spectrum_stat ->SetPoint(7,7.5,0.16213);
    tge_JPsi_forward_spectrum_stat ->SetPoint(8,8.5,0.09443);
    tge_JPsi_forward_spectrum_stat ->SetPoint(9,9.5,0.04853);
    tge_JPsi_forward_spectrum_stat ->SetPoint(10,11,0.02248);

    // statistical error
    tge_JPsi_forward_spectrum_stat ->SetPointError(0,0.5,0.32491);
    tge_JPsi_forward_spectrum_stat ->SetPointError(1,0.5,0.36600);
    tge_JPsi_forward_spectrum_stat ->SetPointError(2,0.5,0.26496);
    tge_JPsi_forward_spectrum_stat ->SetPointError(3,0.5,0.12741);
    tge_JPsi_forward_spectrum_stat ->SetPointError(4,0.5,0.06734);
    tge_JPsi_forward_spectrum_stat ->SetPointError(5,0.5,0.03562);
    tge_JPsi_forward_spectrum_stat ->SetPointError(6,0.5,0.01633);
    tge_JPsi_forward_spectrum_stat ->SetPointError(7,0.5,0.01158);
    tge_JPsi_forward_spectrum_stat ->SetPointError(8,0.5,0.00714);
    tge_JPsi_forward_spectrum_stat ->SetPointError(9,0.5,0.00528);
    tge_JPsi_forward_spectrum_stat ->SetPointError(10,1.0,0.00187);


    // systematic error
    tge_JPsi_forward_spectrum_syst = (TGraphErrors*)tge_JPsi_forward_spectrum_stat->Clone();
    tge_JPsi_forward_spectrum_syst ->SetName("tge_JPsi_forward_spectrum_syst");
    tge_JPsi_forward_spectrum_syst ->SetPointError(0,0.5,0.59610);
    tge_JPsi_forward_spectrum_syst ->SetPointError(1,0.5,0.99704);
    tge_JPsi_forward_spectrum_syst ->SetPointError(2,0.5,0.59763);
    tge_JPsi_forward_spectrum_syst ->SetPointError(3,0.5,0.24519);
    tge_JPsi_forward_spectrum_syst ->SetPointError(4,0.5,0.09296);
    tge_JPsi_forward_spectrum_syst ->SetPointError(5,0.5,0.03838);
    tge_JPsi_forward_spectrum_syst ->SetPointError(6,0.5,0.01744);
    tge_JPsi_forward_spectrum_syst ->SetPointError(7,0.5,0.00721);
    tge_JPsi_forward_spectrum_syst ->SetPointError(8,0.5,0.00465);
    tge_JPsi_forward_spectrum_syst ->SetPointError(9,0.5,0.00236);
    tge_JPsi_forward_spectrum_syst ->SetPointError(10,1.0,0.00115);


    // Normalize
    Double_t integral = 0.0;
    for(Int_t iPoint = 0; iPoint < tge_JPsi_forward_spectrum_stat ->GetN(); iPoint++)
    {
        Double_t pT, y_val;
        Double_t err_X_low  = tge_JPsi_forward_spectrum_stat ->GetErrorXlow(iPoint);
        Double_t err_X_high = tge_JPsi_forward_spectrum_stat ->GetErrorXhigh(iPoint);
        tge_JPsi_forward_spectrum_stat ->GetPoint(iPoint,pT,y_val);

        integral += y_val*(err_X_low + err_X_high);
    }

    if(integral >= 0.0)
    {
        for(Int_t iPoint = 0; iPoint < tge_JPsi_forward_spectrum_stat ->GetN(); iPoint++)
        {
            Double_t pT, y_val;
            tge_JPsi_forward_spectrum_stat ->GetPoint(iPoint,pT,y_val);
            Double_t err_X_low       = tge_JPsi_forward_spectrum_stat ->GetErrorXlow(iPoint);
            Double_t err_X_high      = tge_JPsi_forward_spectrum_stat ->GetErrorXhigh(iPoint);
            Double_t err_Y_low       = tge_JPsi_forward_spectrum_stat ->GetErrorYlow(iPoint);
            Double_t err_Y_high      = tge_JPsi_forward_spectrum_stat ->GetErrorYhigh(iPoint);
            Double_t err_Y_low_syst  = tge_JPsi_forward_spectrum_syst ->GetErrorYlow(iPoint);
            Double_t err_Y_high_syst = tge_JPsi_forward_spectrum_syst ->GetErrorYhigh(iPoint);
            tge_JPsi_forward_spectrum_stat ->SetPoint(iPoint,pT,y_val/integral);
            tge_JPsi_forward_spectrum_syst ->SetPoint(iPoint,pT,y_val/integral);
            tge_JPsi_forward_spectrum_stat ->SetPointError(iPoint,err_X_high,err_Y_high/integral);
            tge_JPsi_forward_spectrum_syst ->SetPointError(iPoint,err_X_high,err_Y_high_syst/integral);
        }
    }
    // End J/Psi
    //----------------------------------------------------------------------




    tgae_dN_dpT_mesons_data[0] = (TGraphAsymmErrors*)vec_tgae_pT_spectra[0][7];    // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
    tgae_dN_dpT_mesons_data[1] = (TGraphAsymmErrors*)vec_tgae_pT_spectra[1][7];    // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
    tgae_dN_dpT_mesons_data[2] = (TGraphAsymmErrors*)vec_tgae_pT_spectra[2][7];    // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
    tgae_dN_dpT_mesons_data[3] = (TGraphAsymmErrors*)tge_phi_dNdpT;    // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
    tgae_dN_dpT_mesons_data[4] = (TGraphAsymmErrors*)tge_Omega_dNdpT[2];    // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
    tgae_dN_dpT_mesons_data[5] = (TGraphAsymmErrors*)tge_D_dNdpT;    // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
    tgae_dN_dpT_mesons_data[6] = (TGraphAsymmErrors*)tge_JPsi_spectra[1][0]; // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
    tgae_dN_dpT_mesons_data[7] = (TGraphAsymmErrors*)tge_D_dNdpT;    // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
}
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
void init_data()
{
    printf("Initialize data \n");
    // https://arxiv.org/pdf/1405.4632.pdf
    inputfile_id   = TFile::Open("./Data/HEPData-ins1297103-v1-root.root"); // https://www.hepdata.net/record/ins1297103
    inputfile_JPsi = TFile::Open("./Data/HEPData-ins1225273-v1-Table_1.root");
    inputfile_D    = TFile::Open("./Data/HEPData-ins1233087-v1-root.root");

    // pi, K+/-, K0s, <K>, p, phi, Lambda, Xi, Omega

    const Int_t arr_centrality_low[7]        = {0,5,10,20,30,40,50};
    const Int_t arr_centrality_high[7]       = {5,10,20,30,40,50,60};
    const Int_t arr_centrality_phi_low[5]    = {10,20,30,40,50};
    const Int_t arr_centrality_phi_high[5]   = {20,30,40,50,60};
    const Int_t arr_centrality_Omega_low[6]  = {5,10,20,30,40,50};
    const Int_t arr_centrality_Omega_high[6] = {10,20,30,40,50,60};
    for(Int_t i = 0; i < 7; i++) // 0..6
    {
        TString label = Form("pi_%d_%d",arr_centrality_low[i],arr_centrality_high[i]);
        arr_labels.push_back(label);
        arr_pid.push_back(0);
    }
    for(Int_t i = 0; i < 7; i++) // 7..13
    {
        TString label = Form("K_%d_%d",arr_centrality_low[i],arr_centrality_high[i]);
        arr_labels.push_back(label);
        arr_pid.push_back(1);
    }
    for(Int_t i = 0; i < 7; i++) // 14..20
    {
        TString label = Form("K0s_%d_%d",arr_centrality_low[i],arr_centrality_high[i]);
        arr_labels.push_back(label);
        arr_pid.push_back(2);
    }
    for(Int_t i = 0; i < 7; i++) // 21..27
    {
        TString label = Form("averK_%d_%d",arr_centrality_low[i],arr_centrality_high[i]);
        arr_labels.push_back(label);
        arr_pid.push_back(3);
    }
    for(Int_t i = 0; i < 7; i++) // 28..34
    {
        TString label = Form("p_%d_%d",arr_centrality_low[i],arr_centrality_high[i]);
        arr_labels.push_back(label);
        arr_pid.push_back(4);
    }
    for(Int_t i = 0; i < 5; i++) // 35..39
    {
        TString label = Form("phi_%d_%d",arr_centrality_phi_low[i],arr_centrality_phi_high[i]);
        arr_labels.push_back(label);
        arr_pid.push_back(5);
    }
    for(Int_t i = 0; i < 7; i++) // 40..46
    {
        TString label = Form("Lambda_%d_%d",arr_centrality_low[i],arr_centrality_high[i]);
        arr_labels.push_back(label);
        arr_pid.push_back(6);
    }
    for(Int_t i = 0; i < 7; i++) // 47..53
    {
        TString label = Form("Xi_%d_%d",arr_centrality_low[i],arr_centrality_high[i]);
        arr_labels.push_back(label);
        arr_pid.push_back(7);
    }
    for(Int_t i = 0; i < 6; i++) // 54..59
    {
        TString label = Form("Omega_%d_%d",arr_centrality_Omega_low[i],arr_centrality_Omega_high[i]);
        arr_labels.push_back(label);
        arr_pid.push_back(8);
    }

    
    for(Int_t i = 0; i < (Int_t)arr_labels.size(); i++)
    {
        vec_graphs.push_back((TGraphAsymmErrors*)inputfile_id->Get(Form("Table %d/Graph1D_y1",i+1)));
        vec_graphs[i]->SetName(arr_labels[i]);
        vec_graphs[i]->SetMarkerColor(arr_color[arr_pid[i]]);
        vec_graphs[i]->SetMarkerStyle(20);
        vec_graphs[i]->SetMarkerSize(0.8);
        vec_graphs[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        vec_graphs[i]->GetYaxis()->SetTitle("v_{2}");
    }


    tg_Upsilon_v2_vs_pT = new TGraphErrors();
    tg_Upsilon_v2_vs_pT ->SetPoint(0,1.88571,0.0129693);
    tg_Upsilon_v2_vs_pT ->SetPoint(1,4.42286,	-0.0109215);
    tg_Upsilon_v2_vs_pT ->SetPoint(2,8.88,	0.00273038);
    tg_Upsilon_v2_vs_pT ->SetPointError(0,0.0,0.040273);
    tg_Upsilon_v2_vs_pT ->SetPointError(1,0.0,0.040273);
    tg_Upsilon_v2_vs_pT ->SetPointError(2,0.0,0.0580205);
    tg_Upsilon_v2_vs_pT ->GetYaxis()->SetRangeUser(-0.1,0.17);
    tg_Upsilon_v2_vs_pT ->SetMarkerStyle(20);
    tg_Upsilon_v2_vs_pT ->SetMarkerSize(0.8);
    tg_Upsilon_v2_vs_pT ->SetMarkerColor(kBlack);

    tg_JPsi_v2_vs_pT = (TGraphAsymmErrors*)inputfile_JPsi->Get("Table 1/Graph1D_y1");
    tg_JPsi_v2_vs_pT ->SetMarkerStyle(20);
    tg_JPsi_v2_vs_pT ->SetMarkerSize(0.8);
    tg_JPsi_v2_vs_pT ->SetMarkerColor(kRed);

    tg_D0_v2_vs_pT = (TGraphAsymmErrors*)inputfile_D->Get("Table 1/Graph1D_y1");
    tg_D0_v2_vs_pT ->SetMarkerStyle(20);
    tg_D0_v2_vs_pT ->SetMarkerSize(0.8);
    tg_D0_v2_vs_pT ->SetMarkerColor(kGray+1);

    //tgae_v2_vs_pT_mesons_data[8]; // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
    tgae_v2_vs_pT_mesons_data[0] = (TGraphAsymmErrors*)vec_graphs[4]->Clone();  // pi
    tgae_v2_vs_pT_mesons_data[1] = (TGraphAsymmErrors*)vec_graphs[18]->Clone(); // K
    tgae_v2_vs_pT_mesons_data[2] = (TGraphAsymmErrors*)vec_graphs[32]->Clone(); // p
    tgae_v2_vs_pT_mesons_data[3] = (TGraphAsymmErrors*)vec_graphs[37]->Clone(); // phi  to be changed
    tgae_v2_vs_pT_mesons_data[4] = (TGraphAsymmErrors*)vec_graphs[57]->Clone(); // Omega  to be changed
    tgae_v2_vs_pT_mesons_data[5] = (TGraphAsymmErrors*)tg_D0_v2_vs_pT->Clone(); // D0
    tgae_v2_vs_pT_mesons_data[6] = (TGraphAsymmErrors*)tg_JPsi_v2_vs_pT->Clone(); // J/Psi
    tgae_v2_vs_pT_mesons_data[7] = (TGraphAsymmErrors*)tg_Upsilon_v2_vs_pT->Clone(); // Upsilon

    for(Int_t i_mass = 0; i_mass < 8; i_mass++)
    {
        tgae_v2_vs_pT_mesons_data_copy[i_mass]  = (TGraphAsymmErrors*)tgae_v2_vs_pT_mesons_data[i_mass] ->Clone();
        tgae_v2_vs_pT_mesons_data_copyB[i_mass] = (TGraphAsymmErrors*)tgae_v2_vs_pT_mesons_data[i_mass] ->Clone();
    }
}
//------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void make_5_60_spectra()
{

    //-----------get points from TGraphs and make 5-60 TGraphs-------------------------------------
    vec_tge_v2_vs_pT_560_pid.resize(3); // pi, K, p

    vector<TString> mini_arr_label = {"pi_5_60","k_5_60","p_5_60"};
    Int_t mini_arr_color[3] = {kBlack,kBlue,kCyan};
    Int_t pid_helper = 0; //because v2 for p are graphs 29..34

    for(Int_t i_pid = 0; i_pid < 3; i_pid++)
    {
        if (i_pid<2) pid_helper = 0;
        if (i_pid==2) pid_helper = 14;

        //-------find max number of v2(pT) points for each particle type------

        Int_t n_arr_max = vec_graphs[1+i_pid*7+pid_helper]->GetN();
        //cout << "n_arr_max:" << n_arr_max << endl;

        for (Int_t i_cent = 1; i_cent < 7; i_cent++)
        {
            //cout << "vec_graphs[i_cent]->GetN():" << vec_graphs[i_cent+i_pid*7+pid_helper]->GetN() << endl;
            //cout << "i_pid:" << i_pid << endl;
            if(n_arr_max > vec_graphs[i_cent+i_pid*7+pid_helper]->GetN())
            {
                n_arr_max = vec_graphs[i_cent+i_pid*7+pid_helper]->GetN();
                //cout << "n_arr_max:" << n_arr_max << endl;
            }
        }

        //-----------------------------------------------------------------

        vec_tge_v2_vs_pT_560_pid[i_pid] = new TGraphErrors(n_arr_max);
        Double_t v2_unnorm_pid[n_arr_max];
        Double_t norm_pid[n_arr_max];
        Double_t x_arr_pid[n_arr_max];
        Double_t v2_pid[n_arr_max];
        Double_t y_err[n_arr_max];

        for(Int_t i_cent = 1; i_cent < 7; i_cent++) // centrality loop
        {

            Int_t n_arr = vec_graphs[i_cent+i_pid*7+pid_helper]->GetN();
            Double_t y_arr_pid,y_pt_arr_pid,v2_y_err;

            //------ make variable arrays 0 before adding centralities

            if(i_cent == 1)
            {
                for(Int_t i_pT = 0; i_pT < n_arr_max; i_pT++)
                {
                    v2_unnorm_pid[i_pT]    = 0;
                    norm_pid[i_pT]         = 0;
                    y_err[i_pT]            = 0;
                }
            }

            //cout << "v2_unnorm_pi[i]:" << v2_unnorm_pi[i] << endl;
            //cout << "norm_pi[i]:" << norm_pi[i] << endl;

            Int_t arr_cent_match[7] = {0,1,2,6,7,8,9};
            for(Int_t i_pT = 0; i_pT < n_arr; i_pT++) // pT loop
            {
                vec_graphs[i_cent+i_pid*7+pid_helper]                     ->GetPoint(i_pT,x_arr_pid[i_pT],y_arr_pid);

                // centrality for pT spectra: 0-5, 5-10, 10-20, 20-40, 40-60, 60-80, 20-30, 30-40, 40-50, 50-60
                // centrality for v2 spectra: 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60
                y_pt_arr_pid        = vec_tgae_pT_spectra[i_pid][arr_cent_match[i_cent]] ->Eval(x_arr_pid[i_pT]);
                v2_unnorm_pid[i_pT] = v2_unnorm_pid[i_pT] + y_arr_pid*y_pt_arr_pid;
                norm_pid[i_pT]      = norm_pid[i_pT] + y_pt_arr_pid;
                v2_y_err            = vec_graphs[i_cent+i_pid*7+pid_helper]  ->GetErrorY(i_pT);
                y_err[i_pT]         = y_err[i_pT] + v2_y_err*y_pt_arr_pid;

                //if(i_pid == 0) printf("i_pid: %d, i_cent: %d, i_pT: %d, pT: %4.2f, v2: %4.3f, dNdpT: %4.6f, v2_unnorm_pid: %4.3f \n",i_pid,i_cent,i_pT,x_arr_pid[i_pT],y_arr_pid,y_pt_arr_pid,v2_unnorm_pid[i_pT]);
            }
        }


        for(Int_t i_pT = 0; i_pT < n_arr_max; i_pT++)
        {
            //cout << "x_arr_pi[i_pT] final:" << x_arr_pi[i_pT] << endl;
            if(norm_pid[i_pT] != 0)
            {
                v2_pid[i_pT] = v2_unnorm_pid[i_pT]/norm_pid[i_pT];
                y_err[i_pT]  = y_err[i_pT]/norm_pid[i_pT];
                vec_tge_v2_vs_pT_560_pid[i_pid] ->SetPoint(i_pT,x_arr_pid[i_pT],v2_pid[i_pT]);
                vec_tge_v2_vs_pT_560_pid[i_pid] ->SetPointError(i_pT,0,y_err[i_pT]);
                //if(i_pid == 0) printf("i_pT: %d, pT: %4.3f, v2_pid: %4.3f \n",i_pT,x_arr_pid[i_pT],v2_pid[i_pT]);
            }

            //cout << "i_pT:" << i_pT << endl;

        }

        //vec_tge_v2_vs_pT_560_pid[i_pid]->SetName("pi_5_60");
        vec_tge_v2_vs_pT_560_pid[i_pid]->SetName(mini_arr_label[i_pid]);
        vec_tge_v2_vs_pT_560_pid[i_pid]->SetMarkerColor(mini_arr_color[i_pid]);
        vec_tge_v2_vs_pT_560_pid[i_pid]->SetMarkerStyle(20);
        vec_tge_v2_vs_pT_560_pid[i_pid]->SetMarkerSize(0.8);
        vec_tge_v2_vs_pT_560_pid[i_pid]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        vec_tge_v2_vs_pT_560_pid[i_pid]->GetYaxis()->SetTitle("v_{2}");
    }

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color

    if((x<0||y<0) && NDC == 1)
    {   // defaults
      x=gPad->GetLeftMargin()*1.15;
      y=(1-gPad->GetTopMargin())*1.04;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextFont(font);
    text->SetTextSize(size);
    if(NDC == 1) text->SetNDC();
    text->SetTextColor(color);
    text->SetTextAngle(angle);
    text->SetTextAlign(align);
    text->Draw();
    return text;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TLine* PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x1_val);
    Zero_line -> SetX2(x2_val);
    Zero_line -> SetY1(y1_val);
    Zero_line -> SetY2(y2_val);
    Zero_line -> SetLineWidth(LineWidth);
    Zero_line -> SetLineStyle(LineStyle);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;

    return Zero_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TCanvas* Draw_1D_histo_and_canvas(TH1D* hist, TString name, Int_t x_size, Int_t y_size,
                              Double_t min_val, Double_t max_val, TString option)
{
    TCanvas* canvas = new TCanvas(name.Data(),name.Data(),10,10,x_size,y_size);
    canvas->SetFillColor(10);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.2);
    canvas->SetRightMargin(0.05);
    canvas->SetLeftMargin(0.2);
    canvas->SetTicks(1,1);
    canvas->SetGrid(0,0);

    hist->SetStats(0);
    hist->SetTitle("");
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetNdivisions(505,'N');
    hist->GetYaxis()->SetNdivisions(505,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    if(min_val != max_val) hist->GetYaxis()->SetRangeUser(min_val,max_val);
    hist->DrawCopy(option.Data());

    return canvas;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TCanvas* Draw_2D_histo_and_canvas(TH2D* hist, TString name, Int_t x_size, Int_t y_size,
                              Double_t min_val, Double_t max_val, TString option)
{
    TCanvas* canvas = new TCanvas(name.Data(),name.Data(),10,10,x_size,y_size);
    canvas->SetFillColor(10);
    canvas->SetTopMargin(0.05);
    canvas->SetBottomMargin(0.2);
    canvas->SetRightMargin(0.22);
    canvas->SetLeftMargin(0.2);
    canvas->SetTicks(1,1);
    canvas->SetGrid(0,0);

    hist->SetStats(0);
    hist->SetTitle("");
    hist->GetXaxis()->SetNdivisions(505,'N');
    hist->GetYaxis()->SetNdivisions(505,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    if(max_val > min_val)
    {
        hist->GetZaxis()->SetRangeUser(min_val,max_val);
    }
    hist->DrawCopy(option.Data());

    return canvas;
}
//----------------------------------------------------------------------------------------




//----------------------------------------------------------------------------------------
Double_t FlowFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t phi, y, v0, v1, v2, v3, v4;
    v0  = par[0];
    v1  = par[1];
    v2  = par[2];
    v3  = par[3];
    v4  = par[4];
    phi = x_val[0];
    y = v0 * (1.0 + 2.0*v1*TMath::Cos(phi) + 2.0*v2*TMath::Cos(2.0*phi)
              + 2.0*v3*TMath::Cos(3.0*phi) + 2.0*v4*TMath::Cos(4.0*phi));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PtFitFunc2_mod(Double_t* x_val, Double_t* par)
{
    Double_t x, y, m0, Temp, Ampl, shift;
    m0    = par[0];
    Temp  = par[1];
    Ampl  = par[2];
    shift = par[3];
    x = x_val[0];
    //y = Ampl*(x-shift)*sqrt((x-shift)*(x-shift)+m0*m0)*TMath::Exp(-(sqrt((x-shift)*(x-shift)+m0*m0)-m0)/Temp);
    //y = Ampl*x*sqrt(x*x+m0*m0)*TMath::Exp(-(sqrt(x*x+m0*m0)-m0)/Temp);
    //y = Ampl*x*sqrt(x*x+m0*m0)*TMath::Exp(-(sqrt(x*x+m0*m0))/Temp);
    //y = Ampl*sqrt(sqrt(x*x+m0*m0))*TMath::Exp(-(sqrt(x*x+m0*m0))/Temp);
    y = Ampl*x*x*TMath::Exp(-sqrt(x*x+m0*m0)/Temp);
    return y;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
Double_t PtFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t pT, y, Temp;
    Temp  = par[0];
    pT = x_val[0];
    y  = pT*TMath::Exp(-pT/Temp);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc(Double_t* x_val, Double_t* par)
{
    // Original function without radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    // pi, K+/-, K0s, <K>, p, phi, Lambda, Xi, Omega
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 200;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/((Double_t)nbins_phi);
    T = par[0];
    rho_0      = par[1];
    rho_a      = par[2];
    s2         = par[3];
    Double_t R = par[4]; // to make it compatible with the R-dependent function

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < (nbins_phi + 1); i++)
    {
        phi   = phi_start + i*delta_phi;
        rho   = rho_0 + rho_a*TMath::Cos(2.0*phi);
        alpha = (pt/T)*TMath::SinH(rho);
        beta  = (mt/T)*TMath::CosH(rho);

        Inte1 += delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        Inte2 += delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc_cout(Double_t* x_val, Double_t* par)
{
    // Original function without radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    // pi, K+/-, K0s, <K>, p, phi, Lambda, Xi, Omega
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 20;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/((Double_t)nbins_phi);
    T = par[0];
    rho_0      = par[1];
    rho_a      = par[2];
    s2         = par[3];
    Double_t R = par[4]; // to make it compatible with the R-dependent function

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < (nbins_phi + 1); i++)
    {
        phi   = phi_start + i*delta_phi;
        rho   = rho_0 + rho_a*TMath::Cos(2.0*phi);
        alpha = (pt/T)*TMath::SinH(rho);
        beta  = (mt/T)*TMath::CosH(rho);

        printf("phi: %4.2f, rho: %4.3f, beta: %4.2f, TMath::BesselK1(beta): %4.5f \n",phi,rho,beta,TMath::BesselK1(beta));

        Inte1 += delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        Inte2 += delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
        printf("v2: %4.3f, nominator: %4.3f, denominator: %4.3f \n",v2,Inte1,Inte2);
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void BlastWaveSimultaneous(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
    Int_t    npfits   = 0;
    Double_t chi2     = 0.0;
    // pi, K+/-, K0s, <K>, p, phi, Lambda, Xi, Omega

    for(Int_t i = 0; i < N_v2_vs_pt_BW; i++) // loop over PIDs
    {
        p[5] = i; // PID
        //Double_t pt_cut = TMath::Sqrt((Mass[i]+mt_m0_cut)*(Mass[i]+mt_m0_cut) - Mass[i]*Mass[i]); // mt-m0 cut
        if(flag_v2_BW_use[i] == 1)
        {
            for(Int_t ix = 0; ix < tgae_v2_stat_BW[i]->GetN(); ix++)
            {
                Double_t x[] = {tgae_v2_stat_BW[i]->GetX()[ix]};
                if(x[0] > arr_pt_low_cut[ix] && x[0] < arr_pt_high_cut[ix])
                {
                    Double_t y      = tgae_v2_stat_BW[i]->GetY()[ix];
                    Double_t ye     = tgae_v2_stat_BW[i]->GetErrorYhigh(ix);
                    ye = 0.01;
                    //cout << "ix = " << ix << ", x = " << x[0] << ", y = " << y << ", ye = " << ye << endl;
                    Double_t bw_val = BlastWaveFitFunc(x,p);
                    //Double_t bw_val = 0.1;
//                    ye += 0.003;
                    Double_t diff   = (y - bw_val)/ye;
                    chi2 += diff*diff;
                    npfits++;
                }
                //else break;
            }
        }
    }

    fval = chi2;
}
//----------------------------------------------------------------------------------------



#if 1
//----------------------------------------------------------------------------------------
void do_minimization()
{
    Double_t par[6]; // T, rho_0, rho_a, s2, R, pid
    // pid: pi, K+/-, K0s, <K>, p, phi, Lambda, Xi, Omega

    par[4] = 1.0;

    Double_t chi2_per_ndf_best = 100000000.0;
    for(Double_t Temp = 90.0; Temp < 150.0; Temp += 10.0)
    {
        printf("Temp: %4.2f \n",Temp);
        par[0] = Temp;
        for(Double_t rho_0 = 1.0; rho_0 < 10.0; rho_0 += 1.0)
        {
            par[1] = rho_0;
            for(Double_t rho_s = 0.05; rho_s < 1.0; rho_s += 0.05)
            {
                par[2] = rho_s;
                for(Double_t s2 = 0.1; s2 < 1.0; s2 += 0.1)
                {
                    Double_t npfits = 0;
                    Double_t chi2   = 0.0;
                    par[3] = s2;
                    for(Int_t ipid = 0; ipid < N_v2_vs_pt_BW; ipid++) // loop over PIDs
                    {
                        par[5] = ipid;
                        if(!flag_v2_BW_use[ipid]) continue;
                        //printf("ipid: %d, npoints: %d \n",ipid,tgae_v2_stat_BW[ipid]->GetN());
                        for(Int_t ix = 0; ix < tgae_v2_stat_BW[ipid]->GetN(); ix++)
                        {
                            Double_t x[] = {tgae_v2_stat_BW[ipid]->GetX()[ix]};
                            //printf("ix: %d, pT: %4.3f \n",ix,x[0]);
                            if(x[0] > arr_pt_low_cut[ix] && x[0] < arr_pt_high_cut[ix])
                            {
                                Double_t y      = tgae_v2_stat_BW[ipid]->GetY()[ix];
                                Double_t ye     = tgae_v2_stat_BW[ipid]->GetErrorYhigh(ix);
                                //ye = 0.01;
                                //cout << "ix = " << ix << ", x = " << x[0] << ", y = " << y << ", ye = " << ye << endl;
                                Double_t bw_val = BlastWaveFitFunc(x,par);
                                //Double_t bw_val = 0.1;
                                //                    ye += 0.003;
                                Double_t diff   = (y - bw_val)/ye;
                                chi2 += diff*diff;
                                npfits += 1.0;
                                //printf("chi2: %4.3f \n",chi2);
                            }
                        }
                    }

                    if((npfits-4.0) > 0.0)
                    {
                        Double_t chi2_per_ndf = chi2/(npfits-4.0);
                        if(chi2_per_ndf < chi2_per_ndf_best)
                        {
                            chi2_per_ndf_best = chi2_per_ndf;
                            printf("chi2_per_ndf_best: %4.3f \n",chi2_per_ndf_best);
                        }
                    }
                }
            }
        }
    }
}
//----------------------------------------------------------------------------------------
#endif



//----------------------------------------------------------------------------------------
TCanvas* Draw_1D_graph_and_canvas(TGraph* hist, TString name, Int_t x_size, Int_t y_size,
                                  Double_t min_val, Double_t max_val, TString option,
                                  Int_t style, Double_t size, Int_t color)
{
    TCanvas* canvas = new TCanvas(name.Data(),name.Data(),10,10,x_size,y_size);
    canvas->SetFillColor(10);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.2);
    canvas->SetRightMargin(0.05);
    canvas->SetLeftMargin(0.2);
    canvas->SetTicks(1,1);
    canvas->SetGrid(0,0);

    hist->SetTitle("");
    hist->SetMarkerSize(size);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(style);
    if(min_val != max_val) hist->GetYaxis()->SetRangeUser(min_val,max_val);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetNdivisions(505,'N');
    hist->GetYaxis()->SetNdivisions(505,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    //for(Int_t i_point = 0; i_point < hist->GetN(); i_point++)
    //{
    //    Double_t x,y;
    //    hist->GetPoint(i_point,x,y);
    //    printf("i_point: %d, x/y: {%4.2f, %4.2f} \n",i_point,x,y);
    //}

    hist->Draw(option.Data());

    return canvas;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t geometric_shape(Double_t R_scale, Double_t phi, Double_t x_fac, Double_t y_fac,
                         Double_t &x_val, Double_t &y_val)
{
    x_val = R_scale*x_fac*TMath::Cos(TMath::DegToRad()*phi);
    y_val = R_scale*y_fac*TMath::Sin(TMath::DegToRad()*phi);
    Double_t radius = TMath::Sqrt(x_val*x_val + y_val*y_val);

    return radius;
}
//----------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------
vector< vector<TVector3> > geometric_shape_perp(Double_t R_scale, Double_t x_fac, Double_t y_fac)
{
    vector< vector<TVector3> > vec_vec_surf_perp;
    vec_vec_surf_perp.resize(2);

    vector<Double_t> x_val;
    vector<Double_t> y_val;
    x_val.resize(3);
    y_val.resize(3);
    TVector3 vec_z_axis(0.0,0.0,1.0);
    for(Int_t phi = 0.0; phi < 360.0; phi += 10.0)
    {
        Int_t counter = 0;
        for(Double_t dphi = -0.1; dphi <= 0.1; dphi += 0.1)
        {
            geometric_shape(R_scale,phi+dphi,x_fac,y_fac,x_val[counter],y_val[counter]);
            counter++;
        }

        TVector3 vec_point_on_surface(x_val[1],y_val[1],0.0);
        vec_vec_surf_perp[0].push_back(vec_point_on_surface);
        TVector3 vec_surf(x_val[2]-x_val[0],y_val[2]-y_val[0],0.0);
        TVector3 vec_surf_perp = vec_surf.Cross(vec_z_axis);
        if(vec_surf_perp.Mag() != 0.0)
        {
            vec_surf_perp *= 1.0/vec_surf_perp.Mag();
        }
        vec_vec_surf_perp[1].push_back(vec_surf_perp);
    }

    return vec_vec_surf_perp;
}
//---------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void get_geometric_shape(TH2D* h2D_geometric_shape_in, Double_t x_fac, Double_t y_fac, Double_t R_scale)
{
    h2D_geometric_shape_in ->Reset();
    Double_t x_pos, y_pos;
    for(Int_t bin_x = 1; bin_x < h2D_geometric_shape_in->GetNbinsX(); bin_x++)
    {
        Double_t x_pos_point = h2D_geometric_shape_in ->GetXaxis()->GetBinCenter(bin_x);
        for(Int_t bin_y = 1; bin_y < h2D_geometric_shape_in->GetNbinsY(); bin_y++)
        {
            Double_t y_pos_point  = h2D_geometric_shape_in ->GetYaxis()->GetBinCenter(bin_y);
            Double_t radius_point = TMath::Sqrt(x_pos_point*x_pos_point + y_pos_point*y_pos_point);
            Double_t phi          = TMath::RadToDeg()*TMath::ATan2(y_pos_point/(y_fac*R_scale),x_pos_point/(x_fac*R_scale));
            Double_t radius       = geometric_shape(R_scale,phi,x_fac,y_fac,x_pos,y_pos);
            if(radius_point <= radius) h2D_geometric_shape_in ->SetBinContent(bin_x,bin_y,1);
        }
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TVector3 get_boost_vector(Double_t R_scale, Double_t x_fac, Double_t y_fac,
                          Double_t x_pos, Double_t y_pos, Double_t rapidity_long,
                          Double_t rho_0, Double_t rho_a
                         )
{
    Double_t R_x   = R_scale*x_fac;
    Double_t R_y   = R_scale*y_fac;
    Double_t r     = TMath::Sqrt(x_pos*x_pos + y_pos*y_pos);
    Double_t phi_s = TMath::ATan2(y_pos,x_pos); // [-Pi,Pi]
    Double_t r_s   = TMath::Sqrt(TMath::Power(r*TMath::Cos(phi_s)/R_x,2.0) + TMath::Power(r*TMath::Sin(phi_s)/R_y,2.0));
    Double_t phi_b = TMath::ATan(TMath::Tan(phi_s)/TMath::Power(R_y/R_x,2.0));
    //Double_t phi_b = TMath::ATan(TMath::Tan(phi_s));
    Double_t phi_b_mod = phi_b;
    //if(fabs(phi_s) > TMath::Pi()/2.0) phi_b_mod += TMath::Pi(); ???
    if(fabs(phi_b_mod) > TMath::Pi()/2.0) phi_b_mod += TMath::Pi();
    if(phi_b_mod > TMath::Pi()) phi_b_mod -= 2.0*TMath::Pi();
    if(phi_b_mod < -TMath::Pi()) phi_b_mod += 2.0*TMath::Pi();


    //printf("pos: {%4.3f, %4.3f}, r: %4.3f, phi_s: %4.3f, r_s: %4.3f, phi_b: %4.3f, phi_b_mod: %4.3f \n",x_pos,y_pos,r,phi_s,r_s,phi_b,phi_b_mod);

    //TVector3 vec_boost(TMath::Cos(phi_b_mod),TMath::Sin(phi_b_mod),0.0);
    TVector3 vec_boost(TMath::Cos(phi_b_mod),TMath::Sin(phi_b_mod),0.0);
    Double_t rapidity = r_s*(rho_0 + rho_a*TMath::Cos(2.0*phi_b_mod));
    //Double_t beta = fabs(TMath::TanH(rapidity))/TMath::CosH(rapidity_long);
    Double_t beta = fabs(TMath::TanH(rapidity));
    vec_boost *= beta;

    return vec_boost;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Optimize_v2(Double_t rho_0_start, Double_t rho_0_stop, Double_t delta_rho_0,
                 Double_t rho_a_start, Double_t rho_a_stop, Double_t delta_rho_a,
                 Double_t R_x_start, Double_t R_x_stop, Double_t delta_R_x,
                 Double_t Temp_start, Double_t Temp_stop, Double_t delta_Temp,
                 Long64_t N_particles, Double_t R_scale,
                 Double_t &rho_0_best, Double_t &rho_a_best, Double_t &R_x_best, Double_t &Temp_best,
                 TProfile* &tp_v2_vs_pT_opt_best_in, TH1D* &h_dNdpT_best_in
                )
{

    Double_t x_pos_point = 0.0;
    Double_t y_pos_point = 0.0;
    Double_t y_fac       = 1.0;

    vector<TProfile*> tp_v2_vs_pT_opt;
    vector<TH1D*>     h_spectra_vs_pT_opt;
    tp_v2_vs_pT_opt.resize(N_masses);
    h_spectra_vs_pT_opt.resize(N_masses);
    for(Int_t i_mass = 0; i_mass < N_masses; i_mass++)
    {
        tp_v2_vs_pT_opt[i_mass]     = new TProfile(Form("tp_v2_vs_pT_opt_%d",i_mass),Form("tp_v2_vs_pT_opt_%d",i_mass),48,0,12.0);
        h_spectra_vs_pT_opt[i_mass]   = new TH1D(Form("h_spectra_vs_pT_opt_%d",i_mass),Form("h_spectra_vs_pT_opt_%d",i_mass),48,0,12.0);
    }

    Double_t chi2_best  = 100000000.0;
    rho_0_best = 0.0;
    rho_a_best = 0.0;
    R_x_best   = 0.0;

    Int_t itt = 0;
    for(Double_t Temp = Temp_start; Temp < Temp_stop; Temp += delta_Temp)
    {
        printf("Temp: %4.3f \n",Temp);
        f_LevyFitFunc ->SetParameter(1,Temp);
        for(Double_t R_x = R_x_start; R_x < R_x_stop; R_x += delta_R_x)
        {
            printf("R_x: %4.3f \n",R_x);
            Double_t x_fac = R_x;
            get_geometric_shape(h2D_geometric_shape,x_fac,1.0,R_scale);
            for(Double_t rho_0 = rho_0_start; rho_0 < rho_0_stop; rho_0 += delta_rho_0)
            {
                printf("  rho_0: %4.3f \n",rho_0);

                for(Double_t rho_a = rho_a_start; rho_a < rho_a_stop; rho_a += delta_rho_a)
                {
                    //for(Int_t i_quark_mass = 0; i_quark_mass < 4; i_quark_mass++)
                    for(Int_t i_quark_mass = 0; i_quark_mass < 3; i_quark_mass++)
                    {
                        Double_t quark_mass = arr_quark_mass_meson[i_quark_mass];
                        f_LevyFitFunc ->SetParameter(0,quark_mass);
                        for(Int_t i_quark = 0; i_quark < N_particles; i_quark++)
                        {
                            h2D_geometric_shape ->GetRandom2(x_pos_point,y_pos_point);

                            // Sample z-rapidity
                            Double_t z_rapidity = (ran.Rndm()-0.5)*8.0; // get bjorken z-rapidity
                            //Double_t z_rapidity = ran.Gaus(0.0,1.2);

                            Double_t beta_z = (TMath::Exp(2.0*z_rapidity) - 1.0)/(TMath::Exp(2.0*z_rapidity) + 1.0); // = TMath::TanH(z_rapidity), z-beta

                            Double_t quark_thermal_phi       = ran.Rndm()*360.0;
                            Double_t quark_thermal_cos_theta = (ran.Rndm()-0.5)*2.0;  // [-1..1]
                            Double_t quark_thermal_theta     = TMath::ACos(quark_thermal_cos_theta); // get theta
                            Double_t quark_thermal_eta       = -TMath::Log(TMath::Tan(quark_thermal_theta/2.0)); // get eta

                            //Double_t quark_pT          = f_LevyFitFunc->GetRandom();
                            Double_t quark_pT                = f_LevyFitFunc->GetRandom()*TMath::Sin(quark_thermal_theta); // sample p, transform to pT

                            TLorentzVector tlv_quark;
                            tlv_quark.SetPtEtaPhiM(quark_pT,quark_thermal_eta,quark_thermal_phi*TMath::DegToRad(),quark_mass); // thermal
                            Double_t pT_thermic = tlv_quark.Pt();

                            TVector3  tv3_boost  = get_boost_vector(R_scale,x_fac,y_fac,x_pos_point,y_pos_point,z_rapidity,rho_0,rho_a);
                            TVector3  tv3_boost_long(0.0,0.0,beta_z); // longitudinal boost vector

                            // boost order is important! -> work in progress
                            //tlv_quark.Boost(tv3_boost);
                            //tlv_quark.Boost(tv3_boost_long);
                            //tlv_quark.Boost(tv3_boost);

#if 1
                            if(ran.Rndm() > 0.1)
                            {
                                //tlv_quark.Boost(tv3_boost);
                                tlv_quark.Boost(tv3_boost_long);
                                tlv_quark.Boost(tv3_boost);
                            }
                            else
                            {
                                tlv_quark.Boost(tv3_boost);
                                tlv_quark.Boost(tv3_boost_long);
                                //tlv_quark.Boost(tv3_boost);
                            }
#endif

                            Double_t pT_lab   = tlv_quark.Pt();
                            Double_t cos_phin = TMath::Cos(2.0*tlv_quark.Phi());
                            Double_t eta      = tlv_quark.Eta();
                            Double_t rapidity = tlv_quark.Rapidity();
                            //if(fabs(eta) > 1.0) continue;
                            if(i_quark_mass < 3  && fabs(eta) > 1.0) continue;
                            if(i_quark_mass < 3  && fabs(rapidity) > 0.5) continue;
                            if(i_quark_mass >= 3 && (rapidity < 2.5 || rapidity > 4.0)) continue;

                            tp_v2_vs_pT_opt[i_quark_mass]    ->Fill(pT_lab,cos_phin);
                            h_spectra_vs_pT_opt[i_quark_mass]  ->Fill(pT_lab);
                        }
                        Double_t integral = h_spectra_vs_pT_opt[i_quark_mass] ->Integral("width");
                        if(integral > 0.0) h_spectra_vs_pT_opt[i_quark_mass] ->Scale(1.0/integral);
                    }

                    // Calculate chi2
                    Double_t chi2 = 0.0;

                    Int_t plot_centrality = 4;


                    //----------------------------------------------------------------------------------------
                    Double_t chi2_per_point_JPsi_spectra = 0.0;
#if 0
                    // J/Psi spectra
                    chi2 = 0.0;
                    Double_t N_used_JPsi_spectra = 0.0;
                    for(Int_t i_point = 0; i_point < tge_JPsi_forward_spectrum_stat->GetN(); i_point++)
                    {
                        Double_t pT_data, dNdpT_data;
                        tge_JPsi_forward_spectrum_stat ->GetPoint(i_point,pT_data,dNdpT_data);
                        if(pT_data > 5.0) break;
                        Double_t dNdpT_err =  tge_JPsi_forward_spectrum_stat->GetErrorYhigh(i_point);
                        Double_t dNdpT_BW  = h_spectra_vs_pT_opt[3] ->GetBinContent(h_spectra_vs_pT_opt[3]->FindBin(pT_data));
                        chi2 += TMath::Power(dNdpT_data - dNdpT_BW,2.0)/TMath::Power(dNdpT_err,2.0);
                        N_used_JPsi_spectra += 1.0;
                        //printf("  i_point: %d, dNdpT_data: %4.3f, dNdpT_BW: %4.3f, sum chi2: %4.3f \n",i_point,dNdpT_data,dNdpT_BW,chi2);
                        //chi2 += fabs(dNdpT_data - dNdpT_BW);
                    }

                    //printf("chi2 spectra: %4.3f \n",chi2);
                    if(N_used_JPsi_spectra > 0.0)
                    {
                        chi2_per_point_JPsi_spectra = chi2/N_used_JPsi_spectra;
                    }
#endif

                    Double_t chi2_per_point_JPsi_v2 = 0.0;
#if 0
                    // J/Psi v2
                    chi2 = 0.0;
                    Double_t N_used_JPsi_v2 = 0.0;
                    for(Int_t i_point = 0; i_point < tg_JPsi_v2_vs_pT->GetN(); i_point++)
                    {
                        Double_t pT_data, v2_data;
                        tg_JPsi_v2_vs_pT ->GetPoint(i_point,pT_data,v2_data);
                        if(pT_data > 5.0) break;
                        Double_t v2_err =  tg_JPsi_v2_vs_pT->GetErrorYhigh(i_point);
                        Double_t v2_BW = tp_v2_vs_pT_opt[3] ->GetBinContent(tp_v2_vs_pT_opt[3]->FindBin(pT_data));
                        chi2 += TMath::Power(v2_data - v2_BW,2.0)/TMath::Power(v2_err,2.0);
                        //printf("  i_point: %d, v2_BW: %4.3f \n",i_point,v2_BW);
                        N_used_JPsi_v2 += 1.0;
                        //chi2 += fabs(v2_data - v2_BW);
                    }

                    //printf("chi2 v2: %4.3f \n",chi2);
                    if(N_used_JPsi_v2 > 0.0)
                    {
                        chi2_per_point_JPsi_v2 = chi2/N_used_JPsi_v2;
                    }

#endif
                    //----------------------------------------------------------------------------------------



                    //----------------------------------------------------------------------------------------
                    Double_t chi2_per_point_pion_spectra = 0.0;
#if 0
                    // Pion spectra
                    chi2 = 0.0;
                    Double_t N_used_pion_spectra = 0.0;
                    for(Int_t i_point = 0; i_point < vec_tgae_pT_spectra[0][3]->GetN(); i_point++)
                    {
                        Double_t pT_data, spectra_data;
                        vec_tgae_pT_spectra[0][3] ->GetPoint(i_point,pT_data,spectra_data);
                        if(pT_data > 1.8) break;
                        Double_t spectra_err = vec_tgae_pT_spectra[0][3] ->GetErrorYhigh(i_point);
                        Double_t spectra_BW = h_spectra_vs_pT_opt[0] ->GetBinContent(h_spectra_vs_pT_opt[0]->FindBin(pT_data));
                        chi2 += TMath::Power(spectra_data - spectra_BW,2.0)/TMath::Power(spectra_err,2.0);
                        //chi2 += fabs(spectra_data - spectra_BW);
                        N_used_pion_spectra += 1.0;
                    }

                    if(N_used_pion_spectra > 0.0)
                    {
                        chi2_per_point_pion_spectra = chi2/N_used_pion_spectra;
                    }
#endif

                    Double_t chi2_per_point_pion_v2 = 0.0;
#if 1
                    // Pion v2
                    chi2 = 0.0;
                    Double_t N_used_pion_v2 = 0.0;
                    for(Int_t i_point = 0; i_point < vec_graphs[plot_centrality]->GetN(); i_point++)
                    {
                        Double_t pT_data, v2_data;
                        vec_graphs[plot_centrality] ->GetPoint(i_point,pT_data,v2_data);
                        if(pT_data > 1.4) break;  // 1.8
                        Double_t v2_err = vec_graphs[plot_centrality] ->GetErrorYhigh(i_point);
                        Double_t v2_BW = tp_v2_vs_pT_opt[0] ->GetBinContent(tp_v2_vs_pT_opt[0]->FindBin(pT_data));
                        //chi2 += TMath::Power(v2_data - v2_BW,2.0)/TMath::Power(v2_err,2.0);
                        chi2 += TMath::Power(fabs(v2_data - v2_BW),2);
                        N_used_pion_v2 += 1.0;
                    }

                    if(N_used_pion_v2 > 0.0)
                    {
                        chi2_per_point_pion_v2 = chi2/N_used_pion_v2;
                    }
#endif
                    //----------------------------------------------------------------------------------------



                    //----------------------------------------------------------------------------------------
                    Double_t chi2_per_point_kaon_spectra = 0.0;
#if 0
                    // Kaon spectra
                    chi2 = 0.0;
                    Double_t N_used_kaon_spectra = 0.0;
                    for(Int_t i_point = 0; i_point < vec_tgae_pT_spectra[1][3]->GetN(); i_point++)
                    {
                        Double_t pT_data, spectra_data;
                        vec_tgae_pT_spectra[1][3] ->GetPoint(i_point,pT_data,spectra_data);
                        if(pT_data > 1.8) break;
                        Double_t spectra_err = vec_tgae_pT_spectra[1][3] ->GetErrorYhigh(i_point);
                        Double_t spectra_BW = h_spectra_vs_pT_opt[1] ->GetBinContent(h_spectra_vs_pT_opt[1]->FindBin(pT_data));
                        chi2 += TMath::Power(spectra_data - spectra_BW,2.0)/TMath::Power(spectra_err,2.0);
                        //chi2 += fabs(spectra_data - spectra_BW);
                        N_used_kaon_spectra += 1.0;
                    }

                    if(N_used_kaon_spectra > 0.0)
                    {
                        chi2_per_point_kaon_spectra = chi2/N_used_kaon_spectra;
                    }
#endif

                    Double_t chi2_per_point_kaon_v2 = 0.0;
#if 1
                    // Kaon v2
                    chi2 = 0.0;
                    Double_t N_used_kaon_v2 = 0.0;
                    for(Int_t i_point = 0; i_point < vec_graphs[plot_centrality+14]->GetN(); i_point++)
                    {
                        Double_t pT_data, v2_data;
                        vec_graphs[plot_centrality+14] ->GetPoint(i_point,pT_data,v2_data);
                        if(pT_data > 1.8) break;
                        Double_t v2_err = vec_graphs[plot_centrality+14] ->GetErrorYhigh(i_point);
                        Double_t v2_BW = tp_v2_vs_pT_opt[1] ->GetBinContent(tp_v2_vs_pT_opt[1]->FindBin(pT_data));
                        //chi2 += TMath::Power(v2_data - v2_BW,2.0)/TMath::Power(v2_err,2.0);
                        chi2 += TMath::Power(fabs(v2_data - v2_BW),2);
                        N_used_kaon_v2 += 1.0;
                    }

                    if(N_used_kaon_v2 > 0.0)
                    {
                        chi2_per_point_kaon_v2 = chi2/N_used_kaon_v2;
                    }
#endif
                    //----------------------------------------------------------------------------------------



                    //----------------------------------------------------------------------------------------
                    Double_t chi2_per_point_proton_spectra = 0.0;
#if 0
                    // Proton spectra
                    chi2 = 0.0;
                    Double_t N_used_proton_spectra = 0.0;
                    for(Int_t i_point = 0; i_point < vec_tgae_pT_spectra[2][3]->GetN(); i_point++)
                    {
                        Double_t pT_data, spectra_data;
                        vec_tgae_pT_spectra[2][3] ->GetPoint(i_point,pT_data,spectra_data);
                        if(pT_data > 2.4) break;
                        Double_t spectra_err = vec_tgae_pT_spectra[2][3] ->GetErrorYhigh(i_point);
                        Double_t spectra_BW = h_spectra_vs_pT_opt[2] ->GetBinContent(h_spectra_vs_pT_opt[2]->FindBin(pT_data));
                        chi2 += TMath::Power(spectra_data - spectra_BW,2.0)/TMath::Power(spectra_err,2.0);
                        //chi2 += fabs(spectra_data - spectra_BW);
                        N_used_proton_spectra += 1.0;
                    }

                    if(N_used_proton_spectra > 0.0)
                    {
                        chi2_per_point_proton_spectra = chi2/N_used_proton_spectra;
                    }
#endif

                    Double_t chi2_per_point_proton_v2 = 0.0;
#if 1
                    // proton v2
                    chi2 = 0.0;
                    Double_t N_used_proton_v2 = 0.0;
                    for(Int_t i_point = 0; i_point < vec_graphs[plot_centrality+28]->GetN(); i_point++)
                    {
                        Double_t pT_data, v2_data;
                        vec_graphs[plot_centrality+28] ->GetPoint(i_point,pT_data,v2_data);
                        if(pT_data > 2.4) break;
                        Double_t v2_err = vec_graphs[plot_centrality+28] ->GetErrorYhigh(i_point);
                        Double_t v2_BW = tp_v2_vs_pT_opt[2] ->GetBinContent(tp_v2_vs_pT_opt[2]->FindBin(pT_data));
                        //chi2 += TMath::Power(v2_data - v2_BW,2.0)/TMath::Power(v2_err,2.0);
                        chi2 += TMath::Power(fabs(v2_data - v2_BW),2);
                        N_used_proton_v2 += 1.0;
                    }

                    if(N_used_proton_v2 > 0.0)
                    {
                        chi2_per_point_proton_v2 = chi2/N_used_proton_v2;
                    }
#endif
                    //----------------------------------------------------------------------------------------



                    Double_t chi2_per_point_total = chi2_per_point_JPsi_spectra + chi2_per_point_JPsi_v2
                        + chi2_per_point_pion_spectra + chi2_per_point_pion_v2
                        + chi2_per_point_kaon_spectra + chi2_per_point_kaon_v2
                        + chi2_per_point_proton_spectra + chi2_per_point_proton_v2;

                    if(chi2_per_point_total  < chi2_best)
                    {
                        chi2_best  = chi2_per_point_total;
                        rho_0_best = rho_0;
                        rho_a_best = rho_a;
                        R_x_best   = R_x;
                        Temp_best  = Temp;
                        if(tp_v2_vs_pT_opt_best_in) tp_v2_vs_pT_opt_best_in ->Delete();
                        tp_v2_vs_pT_opt_best_in = (TProfile*)tp_v2_vs_pT_opt[0]->Clone("tp_v2_vs_pT_opt_best_in");

                        if(h_dNdpT_best_in) h_dNdpT_best_in ->Delete();
                        h_dNdpT_best_in = (TH1D*)h_spectra_vs_pT_opt[3]->Clone("h_dNdpT_best_in");

                        printf("itt: %d, chi2_best: %4.5f, Temp: %4.4f, rho_0: %4.3f, rho_a: %4.3f, R_x: %4.3f \n",itt,chi2_best,Temp_best,rho_0_best,rho_a_best,R_x_best);
                    }

                    for(Int_t i_mass = 0; i_mass < N_masses; i_mass++)
                    {
                        tp_v2_vs_pT_opt[i_mass] ->Reset();
                    }


                    itt++;
                }
            }
        }
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t get_s2(Double_t R_x, Double_t R_y)
{
    // Calculates s2 from Rx and Rx, based on equation (6) in
    // https://arxiv.org/pdf/nucl-th/0312024.pdf
    // WARNING: Only corrrect if rho2 is 0.

    Double_t s2 = 0.5*(TMath::Power(R_y/R_x,2.0) - 1.0)/(TMath::Power(R_y/R_x,2.0) + 1.0);

    return s2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_Circle_2D_new(Float_t radius_in = 1.0, Float_t radius_out = 2.0,const Int_t n_radii = 1,
                        const Int_t n_delta_phi = 2, Float_t color = 2, Int_t line_style = 1, Int_t line_width = 1,
                        Float_t x_offset = 0.0, Float_t y_offset = 0.0, Int_t fill_color = 1,
                        Double_t transparency_line = 0.0, Double_t transparency_fill = 0.0
                       )
{
    Float_t z = 0.0;
    const Int_t n_points = 50;
    TPolyLine   *tp_Circles[n_radii];
    TPolyLine   *tp_Radial[n_delta_phi];
    Float_t radius_table[n_radii];
    Float_t delta_radius;
    if(n_radii > 1) {delta_radius = (radius_out-radius_in)/((Float_t)(n_radii-1));}
    else{delta_radius = 0.0;}
    Float_t delta_phi    = 2.0*TMath::Pi()/((Float_t)n_delta_phi);
    Float_t z_tpc_val    = z;

    for(Int_t r = 0; r < n_radii; r++)
    {
        radius_table[r] = radius_in + r*delta_radius;
        tp_Circles[r] = new TPolyLine();
        Float_t radius   = radius_table[r];
        for(Int_t t = 0; t < n_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            //cout << ", x: " << x_tpc_val << ", y: " << y_tpc_val << endl;
            tp_Circles[r]->SetNextPoint(x_tpc_val,y_tpc_val);
        }

        tp_Circles[r]->SetLineStyle(line_style);
        tp_Circles[r]->SetLineColorAlpha(color,transparency_line); // 28
        tp_Circles[r]->SetLineWidth(line_width);
        tp_Circles[r]->SetFillColorAlpha(fill_color,transparency_fill);
        //tp_Circles[r]->DrawClone("oglf");
        tp_Circles[r]->Draw("f");
        tp_Circles[r]->Draw("ogl");
    }

    for(Int_t r = 0; r < n_delta_phi; r++)
    {
        tp_Radial[r] = new TPolyLine();
        Float_t phi_val = r*delta_phi;
        for(Int_t t = 0; t < 2; t++)
        {
            Float_t radius;
            if(t == 0) {radius = radius_table[0];}
            else {radius = radius_table[n_radii-1];}
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            tp_Radial[r]->SetNextPoint(x_tpc_val,y_tpc_val);
        }
        tp_Radial[r]->SetLineStyle(0);
        tp_Radial[r]->SetLineColor(color); // 28
        tp_Radial[r]->SetLineWidth(1);
        tp_Radial[r]->Draw("ogl");
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Init_v2_Mathematica()
{
    // From Klaus Reygers -> Mathematica
    printf("Init_v2_Mathematica \n");

    vec_tg_v2_vs_pT_Mathematica.resize(5); // pi, K, p, J/Psi, Upsilon
    vec_tg_dNdpT_vs_pT_Mathematica.resize(5); // pi, K, p, J/Psi, Upsilon

    vector<std::fstream> myfile_v2;
    myfile_v2.resize(5); // pi, K, p, J/Psi, Upsilon
    myfile_v2[0].open((char*)"./Data/v2_pion_klaus.csv");
    myfile_v2[1].open((char*)"./Data/v2_kaon_klaus.csv");
    myfile_v2[2].open((char*)"./Data/v2_proton_klaus.csv");
    myfile_v2[3].open((char*)"./Data/v2_jpsi_klaus.csv");
    myfile_v2[4].open((char*)"./Data/v2_upsilon_klaus.csv");

    vector<std::fstream> myfile_dNdpT;
    myfile_dNdpT.resize(5); // pi, K, p, J/Psi, Upsilon
    myfile_dNdpT[0].open((char*)"./Data/inv_yield_pion_klaus.csv");
    myfile_dNdpT[1].open((char*)"./Data/inv_yield_kaon_klaus.csv");
    myfile_dNdpT[2].open((char*)"./Data/inv_yield_proton_klaus.csv");
    myfile_dNdpT[3].open((char*)"./Data/inv_yield_jpsi_klaus.csv");
    myfile_dNdpT[4].open((char*)"./Data/inv_yield_upsilon_klaus.csv");

    for(Int_t i_file = 0; i_file < (Int_t)vec_tg_v2_vs_pT_Mathematica.size(); i_file++)
    {
        vec_tg_v2_vs_pT_Mathematica[i_file]    = new TGraph();
        vec_tg_dNdpT_vs_pT_Mathematica[i_file] = new TGraph();
        Int_t N_points = 0;
        //printf("Open i_file: %d \n",i_file);
        while(!myfile_v2[i_file].eof())
        {
            //cout << "file: " << i_file << ", line: " << N_points << endl;
            //if(line_counter > 10) break;
            std::string str;
            std::getline(myfile_v2[i_file], str);
            if(str == "") continue;

            std::size_t found = str.find(",");
            std::string sub_str;
            sub_str = str.substr(0,found);
            Double_t pT = std::stod(sub_str);
            str.replace(0,found+1,"");

            Double_t v2 = std::stod(str);

            //printf("i_file: %d, pT: %4.3f, v2: %4.3f \n",i_file,pT,v2);

            vec_tg_v2_vs_pT_Mathematica[i_file] ->SetPoint(N_points,pT,v2);
            N_points++;
        }

        Double_t integral = 0.0;
        Double_t pT_A = 0.0;
        Double_t pT_B = 0.0;
        N_points = 0;
        while(!myfile_dNdpT[i_file].eof())
        {
            //cout << "file: " << i_file << ", line: " << N_points << endl;
            //if(line_counter > 10) break;
            std::string str;
            std::getline(myfile_dNdpT[i_file], str);
            if(str == "") continue;

            std::size_t found = str.find(",");
            std::string sub_str;
            sub_str = str.substr(0,found);
            Double_t pT = std::stod(sub_str);
            str.replace(0,found+1,"");

            Double_t yield = std::stod(str);
            yield *= pT;

            //printf("i_file: %d, pT: %4.3f, yield: %4.3f \n",i_file,pT,yield);

            integral += yield;
            vec_tg_dNdpT_vs_pT_Mathematica[i_file] ->SetPoint(N_points,pT,yield);
            N_points++;
            if(N_points == 1) pT_A = pT;
            if(N_points == 2) pT_B = pT;
            //printf("N_points: %d, pT_A: %4.3f, pT_B: %4.3f \n",N_points,pT_A,pT_B);
        }
        Double_t bin_width = fabs(pT_B - pT_A);
        printf("integral: %4.5f, bin_width: %4.3f \n",integral,bin_width);
        integral *= bin_width;
        if(integral > 0.0)
        {
            for(Int_t i_point = 0; i_point < vec_tg_dNdpT_vs_pT_Mathematica[i_file]->GetN(); i_point++)
            {
                Double_t x, y;
                vec_tg_dNdpT_vs_pT_Mathematica[i_file] ->GetPoint(i_point,x,y);
                vec_tg_dNdpT_vs_pT_Mathematica[i_file] ->SetPoint(i_point,x,y/integral);
            }
        }
    }
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
void Init_density()
{
    // From Klaus Reygers -> Mathematica, Glauber overlap
    h2D_density_Glauber = new TH2D("h2D_density_Glauber","h2D_density_Glauber",201,-10.0,10.0,201,-10.0,10.0);

    std::fstream myfile;
    myfile.open((char*)"./Data/tab_profile_PbPb_b_eq_9fm.csv");
    cout << "Density file opened" << endl;
    Int_t line_counter = 0;
    Double_t max_density = 0.0;
    Double_t xy_range[2][2] = {0.0};
    while(!myfile.eof())
    {
        //if(line_counter > 10) break;
        std::string str;
        std::getline(myfile, str);
        if(str == "") continue;

        std::size_t found = str.find(",");
        std::string sub_str;
        sub_str = str.substr(0,found);
        Double_t x_val = std::stod(sub_str);
        str.replace(0,found+1,"");

        found = str.find(",");
        sub_str = str.substr(0,found);
        Double_t y_val = std::stod(sub_str);
        str.replace(0,found+1,"");

        Double_t density = std::stod(str);
        if(density > max_density) max_density = density;

        if(line_counter == 0)
        {
            xy_range[0][0] = x_val; // x_min
            xy_range[0][1] = x_val; // x_max
            xy_range[1][0] = y_val; // y_min
            xy_range[1][1] = y_val; // y_max
        }
        else
        {
            if(x_val < xy_range[0][0]) xy_range[0][0] = x_val;
            if(x_val > xy_range[0][1]) xy_range[0][1] = x_val;
            if(y_val < xy_range[1][0]) xy_range[1][0] = y_val;
            if(y_val > xy_range[1][1]) xy_range[1][1] = y_val;
        }

        //density = TMath::Power(density,5);
        h2D_density_Glauber ->Fill(x_val,y_val,density);

        //printf("pos: {%4.3f, %4.3f, density: %4.16f} \n",x_val,y_val,density);

        line_counter++;
    }

    printf("max_density: %4.6f, x-range: {%4.3f, %4.3f}, y-range: {%4.3f, %4.3f} \n",max_density,xy_range[0][0],xy_range[0][1],xy_range[1][0],xy_range[1][1]);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_hist_line(TH1D* hist, Double_t x_start, Double_t x_stop,
                    Double_t y_start, Double_t y_stop,
                    Double_t color, Int_t width, Int_t style, Double_t trans,
                    TString option)
{
    TPolyLine* line = new TPolyLine();
    for(Int_t i_bin = 1; i_bin <= hist->GetNbinsX(); i_bin++)
    {
        Double_t bin_cont = hist ->GetBinContent(i_bin);
        Double_t bin_cent = hist ->GetBinCenter(i_bin);
        if(bin_cent < x_start || bin_cent > x_stop) continue;
        if(bin_cont < y_start) continue;
        if(bin_cont > y_stop) break;
        line ->SetNextPoint(bin_cent,bin_cont);
    }

    line ->SetLineStyle(style);
    line ->SetLineColorAlpha(color,trans);
    line ->SetLineWidth(width);
    line ->Draw(option.Data());
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void plot_spectra()
{
    // pi, K, p
    // 0-5, 5-10, 10-20, 20-40, 40-60, 60-80, 20-30, 30-40, 40-50

    TCanvas* can_dNdpT_vs_pT = new TCanvas("can_dNdpT_vs_pT","can_dNdpT_vs_pT",10,10,1250,700);
    can_dNdpT_vs_pT ->Divide(3,2);
    for(Int_t iPad = 1; iPad <= 6; iPad++)
    {
        can_dNdpT_vs_pT ->cd(iPad)->SetLogy(1);
        can_dNdpT_vs_pT ->cd(iPad);
    }

    for(Int_t iPid = 0; iPid < 5; iPid++)
    {
        Double_t integral = h_dN_dpT_mesons[iPid] ->Integral("width");
        if(integral <= 0.0) continue;
        h_dN_dpT_mesons[iPid] ->Scale(1.0/integral);
    }


    for(Int_t i_mass = 0; i_mass < 5; i_mass++)
    {
        vec_tg_dNdpT_vs_pT_Mathematica[i_mass] ->SetMarkerColor(arr_color_mass[i_mass]);
        vec_tg_dNdpT_vs_pT_Mathematica[i_mass] ->SetMarkerStyle(28); // 24
        vec_tg_dNdpT_vs_pT_Mathematica[i_mass] ->SetLineColor(arr_color_mass[i_mass]);
        vec_tg_dNdpT_vs_pT_Mathematica[i_mass] ->SetLineStyle(1);
        vec_tg_dNdpT_vs_pT_Mathematica[i_mass] ->SetLineWidth(2);
    }


    // Pions
    can_dNdpT_vs_pT ->cd(1);
    vec_tgae_pT_spectra[0][10] ->Draw("AP");
    h_dN_dpT_mesons[0] ->DrawCopy("same");
    vec_tg_dNdpT_vs_pT_Mathematica[0] ->Draw("same P");


    for(Int_t i_par = 0; i_par < 4; i_par++)
    {
        f_FitBessel->ReleaseParameter(i_par);
        f_FitBessel->SetParError(i_par,0.0);
        f_FitBessel->SetParameter(i_par,0.0);
    }
    f_FitBessel ->FixParameter(0,arr_quark_mass_meson[0]);
    f_FitBessel ->FixParameter(1,0.14);
    f_FitBessel ->FixParameter(3,0.0);
    f_FitBessel ->SetRange(0.0,5.0);
    vec_tg_dNdpT_vs_pT_Mathematica[0] ->Fit("f_FitBessel","WMN","",0.0,5.0);

    f_FitBessel ->SetLineColor(kRed);
    f_FitBessel ->SetLineStyle(1);
    f_FitBessel ->SetLineWidth(3);
    f_FitBessel ->DrawCopy("same");


    // Kaons
    can_dNdpT_vs_pT ->cd(2);
    vec_tgae_pT_spectra[1][10] ->Draw("AP");
    h_dN_dpT_mesons[1] ->DrawCopy("same");
    vec_tg_dNdpT_vs_pT_Mathematica[1] ->Draw("same P");

    for(Int_t i_par = 0; i_par < 4; i_par++)
    {
        f_FitBessel->ReleaseParameter(i_par);
        f_FitBessel->SetParError(i_par,0.0);
        f_FitBessel->SetParameter(i_par,0.0);
    }
    f_FitBessel ->FixParameter(0,arr_quark_mass_meson[1]);
    f_FitBessel ->FixParameter(1,0.14);
    f_FitBessel ->FixParameter(3,0.0);
    f_FitBessel ->SetRange(0.0,5.0);
    vec_tg_dNdpT_vs_pT_Mathematica[1] ->Fit("f_FitBessel","WMN","",0.0,5.0);

    f_FitBessel ->SetLineColor(kRed);
    f_FitBessel ->SetLineStyle(1);
    f_FitBessel ->SetLineWidth(3);
    f_FitBessel ->DrawCopy("same");

    // Protons
    can_dNdpT_vs_pT ->cd(3);
    vec_tgae_pT_spectra[2][10] ->Draw("AP");
    h_dN_dpT_mesons[2] ->DrawCopy("same");
    vec_tg_dNdpT_vs_pT_Mathematica[2] ->Draw("same P");

    for(Int_t i_par = 0; i_par < 4; i_par++)
    {
        f_FitBessel->ReleaseParameter(i_par);
        f_FitBessel->SetParError(i_par,0.0);
        f_FitBessel->SetParameter(i_par,0.0);
    }
    f_FitBessel ->FixParameter(0,arr_quark_mass_meson[2]);
    f_FitBessel ->FixParameter(1,0.14);
    f_FitBessel ->FixParameter(3,0.0);
    f_FitBessel ->SetRange(0.0,5.0);
    vec_tg_dNdpT_vs_pT_Mathematica[2] ->Fit("f_FitBessel","WMN","",0.0,5.0);

    f_FitBessel ->SetLineColor(kRed);
    f_FitBessel ->SetLineStyle(1);
    f_FitBessel ->SetLineWidth(3);
    f_FitBessel ->DrawCopy("same");

    // J/Psi
    can_dNdpT_vs_pT ->cd(4);
    tge_JPsi_spectra[1][0]->Draw("AP"); // 0-20%, 20-40%, 40-90%
    h_dN_dpT_mesons[3] ->DrawCopy("same");
    tge_JPsi_forward_spectrum_stat ->SetLineColor(kRed);
    tge_JPsi_forward_spectrum_stat ->Draw("same");
    vec_tg_dNdpT_vs_pT_Mathematica[3] ->Draw("same P");
    if(h_dNdpT_best)
    {
        h_dNdpT_best ->SetLineColor(kGreen);
        h_dNdpT_best ->DrawCopy("same h");
    }

    // Upsilons
    can_dNdpT_vs_pT ->cd(5);
    vec_tgae_pT_spectra[2][3] ->Draw("AP");
    h_dN_dpT_mesons[4] ->DrawCopy("same");
    vec_tg_dNdpT_vs_pT_Mathematica[4] ->Draw("same P");

    // D0
    can_dNdpT_vs_pT ->cd(6);
    tge_D_dNdpT ->Draw("AP");


    for(Int_t iPad = 1; iPad <= 5; iPad++)
    {
        can_dNdpT_vs_pT ->cd(iPad);
        plotTopLegend((char*)label_pid_spectra[iPad-1].Data(),0.75,0.83,0.06,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        plotTopLegend((char*)"|y|<0.5",0.75,0.77,0.06,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    }
}
//----------------------------------------------------------------------------------------




