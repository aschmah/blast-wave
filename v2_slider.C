
#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"
#include "TGLayout.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGTextEntry.h"
#include "TGTripleSlider.h"
#include <TGSlider.h>

//#include "./feeddown/feeddown.h"
#include "./functions_BW.h"

//R__ADD_LIBRARY_PATH(./feeddown/)
//R__LOAD_LIBRARY(libfeeddown.so)

ClassImp(Tblastwave_yield_and_v2)

enum ETestCommandIdentifiers {
    HId1,
    HId1a,
    HId2,
    HId3,
    HCId1,
    HCId2,
    HCId3,
    HSId1
};
class TBlastWaveGUI : public TGMainFrame {
private:
    TRootEmbeddedCanvas *fCanvas;
    TRootEmbeddedCanvas *fCanvasB;
    TGLayoutHints       *fLcan;
    TGLayoutHints       *fLcanB;
    TF1                 *fFitFcn;
    TGHorizontalFrame   *fHframe0, *fHframe1, *fHframe1a, *fHframe2;
    TGLayoutHints       *fBly, *fBfly1, *fBfly2, *fBfly3;
    TGTripleHSlider     *fHslider1;
    TGHSlider           *hslider;
    TGTextEntry         *fTeh1, *fTeh2, *fTeh3;
    TGTextBuffer        *fTbh1, *fTbh1a, *fTbh2, *fTbh3;
    TGGroupFrame        *fGroupFrames[7];
    TGGroupFrame        *pT_Range_Group[2][N_masses_all];
    TGGroupFrame        *GroupSlider[5];
    TGGroupFrame        *fGroupFrame_PID_fit[2];
    TGGroupFrame        *fGroupFrame_PID_plot[2];
    TGGroupFrame        *fGroupFrames_particles[4][N_masses_all];

    TProfile* tp_v2_vs_pT_mesons[N_masses][9][9][9][9][9] = {NULL}; // [i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]
    TH1F*     h_dN_dpT_mesons[N_masses][9][9][9][9][9]    = {NULL}; // [i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]
    TList* arr_lists[2][9][9];
    vector<TGHorizontalFrame*>  vec_Hframe;
    vector<TGHSlider*>          vec_slider;
    vector<TGLayoutHints*>      vec_LayoutHints;
    vector<TGTextEntry*>        vec_TextEntry;
    vector<TGTextBuffer*>       vec_TextBuffer;
    vector<TGLabel*>            vec_TGLabel;

    TGMainFrame* FrameA;
    TGMainFrame* FrameB;
    TGLayoutHints* fLcanC;

    TGMainFrame* FrameD;
    TGLayoutHints* fLcanD;
    TGHorizontalFrame *hframeD1;
    TGHorizontalFrame *hframeD2;
    TGHorizontalFrame *hframeD2a;
    TGHorizontalFrame *hframeD2b;
    //TGVerticalFrame *hframeD3;
    TGHorizontalFrame *hframeD3;
    TGHorizontalFrame *hframeD4;
    TGHorizontalFrame *hframeD5;
    TGHorizontalFrame *hframeD5a;
    TGHorizontalFrame *hframeD6;
    TGHorizontalFrame *Hframe_pT_limits;
    TGVerticalFrame *hVframeD5a[4];
    TGVerticalFrame *hVframeD3;
    TGVerticalFrame *hVframeD4;
    TGVerticalFrame *Vframe_pT_limits[2];
    TGVerticalFrame *vframeD1;
    TGCompositeFrame *cframe2;
    TGCompositeFrame *fCompositeFrame_pid[8];
    //TGCompositeFrame *fCompositeFrame_pid_plot[2][2];
    TGTextButton *ButtonD1a;
    TGTextButton *Button_exit;
    TGTextButton *Button_save;
    TGHProgressBar* fHProg1 = NULL;
    TGHProgressBar* fHProg2 = NULL;
    TGVProgressBar* fVProg1;
    TGLayoutHints* fHint2;
    TGTextButton      *fGO;
    TGLayoutHints* fHint3;
    TGLayoutHints* LHintsD4a, *LHintsD4a2;
    TGNumberEntry* NEntryD3a;
    TGLabel*       LabelD4a;
    TGNumberEntry* NEntryD3b;
    TGLabel*       LabelD4b;
    TGTransientFrame* frame_TGTransient;
    TGTransientFrame* frame_TGTransientB;
    TGTextButton      *Button_take_params_MC_to_ana;
    TGTextButton      *Button_take_params_Set_to_ana;
    TGComboBox        *fCombo, *ComboEnergy, *ComboCentrality;

    TGCheckButton* fCheckBox_sel[3];
    TGCheckButton* fCheckBox_pid[N_masses_all];
    TGCheckButton* fCheckBox_pid_fit_dNdpt[N_masses_all];
    TGCheckButton* fCheckBox_pid_plot[N_masses_all];
    TGCheckButton* fCheckBox_pid_plot_dNdpt[N_masses_all];
    TGCheckButton* fCheckBox_pid_set[4][N_masses_all];
    TGCheckButton* fCheckBox_v2_dNdpT[2];
    TGCheckButton* fCheckBoxFeedDown;
    TGLayoutHints* fLCheckBox;
    TGLayoutHints* fLCheckBoxB;
    TGLayoutHints* fLGroupFrames_particles;
    TString label_checkbox[3] = {"Plot MC","Plot Ana","Plot chi^2"};

    TGHorizontalFrame* arr_HFrame_NEntry_limits[N_masses];
    TGVerticalFrame* arr_VFrame_NEntry_limits[2];
    TGLabel*           arr_Label_NEntry_limits[2][N_masses_all];
    TGLabel*           arr_Label_NEntry_limits_A[2][N_masses_all];
    TGNumberEntry*     arr_NEntry_limits[2][N_masses_all];
    TGNumberEntry*     arr_NEntry_limits_A[2][N_masses_all];
    TGNumberEntry*     NEntry_set_limits;;
    TGNumberEntry*     arr_NEntry_ana_params[4];
    TGLabel*           arr_Label_NEntry_ana_params[4];
    TGLayoutHints* arr_fL1[2];
    Double_t min_max_pT_range_pid[2][N_masses_all];
    TGraph* tg_v2_BW_ana_pid[N_masses];
    TGraph* tg_dNdpT_BW_ana_pid[N_masses];

    TGraph* tg_v2_BW_ana_pid_min[N_masses_all];
    TGraph* tg_dNdpT_BW_ana_pid_min[N_masses_all];

    TGraph* tg_v2_BW_ana_pid_plot[N_masses];
    TGraph* tg_v2_BW_ana_pid_plot_range[N_masses];
    TGraph* tg_dNdpT_BW_ana_pid_plot[N_masses];
    TGraph* tg_dNdpT_BW_ana_pid_plot_range[N_masses];

    TGraph *gr2;
    TGraph *gr1;


    Double_t var_test = 5.2;

    TGTextButton *Button_minimize;
    TGTextButton *Button_stop_minimize;

    TGTextButton *Button_minimize_ana;
    TGTextButton *Button_stop_minimize_ana;
    TGTextButton *Button_draw_ellipse_ana;

    TGTextButton *Button_make_plot_v2;
    TGTextButton *Button_make_plot_dNdpT;
    TGTextButton *Button_plot_data;
    TGTextButton *Button_write_params;
    TGTextButton *Set_energy_centrality;
    TGTextEntry  *Add_output_file_name;

    TFile* inputfiles[2];
    TH1D* h_dummy;
    TH1D* h_dummy_plot_v2[3] = {NULL};
    TH1D* h_dummy_dNdpT;
    TH1D* h_dummy_plot_data;
    TH1D* h_fit_params;
    TH1D* h_min_val_pT;
    TH1D* h_max_val_pT;
    vector<Double_t> vec_min_val_pT;
    vector<Double_t> vec_max_val_pT;
    TLegend* leg_v2_vs_pT_A = NULL;
    TLegend* leg_v2_vs_pT_B = NULL;
    TLegend* leg_dNdpT_vs_pT = NULL;
    TLegend* leg_v2_vs_pT_C = NULL;
    TLegend* legend_v2_plot_data = NULL;
    TLegend* legend_dNdpt_plot_data = NULL;
    TGraph* tg_leg = NULL;
    Int_t flag_stop_minimize = 0;
    Int_t flag_minimization_ana = 0;
    Double_t integration_range_pid[N_masses][2] = {0.0};
    Double_t individual_chi2_ana[N_masses_all][2] = {0.0};
    Double_t best_individual_chi2_per_point_ana[N_masses][2] = {0.0};
    Double_t N_individual_chi2_ana[N_masses_all][2] = {0.0};

    Double_t chi2_min;
    Double_t chi2_ndf_dNdpT_min;
    Double_t chi2_final_min = 0.0;
    Double_t chi2_final_min_ana = 999999.0;
    Double_t chi2_final_min_piKp;
    TString HistName;
    char NoP[50];

    TLatex* TLatex_chi2_MC = NULL;
    TLatex* TLatex_chi2_ana = NULL;
    TLatex* TLatex_legend_system[4] = {NULL,NULL,NULL,NULL};
    TLatex* TLatex_legend_dNdpT[N_masses] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TLatex* TLatex_legend_v2_plot[N_masses] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TLatex* TLatex_legend_dNdpT_plot[N_masses] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    TLine* TL_line_base = NULL;
    TLine* TL_line_base_plot[3] = {NULL};

    Double_t T_BW_fit_ana        = -1.0;
    Double_t Rho0_BW_fit_ana     = -1.0;
    Double_t Rho2_BW_fit_ana     = -1.0;
    Double_t RxOverRy_BW_fit_ana = -1.0;

    TArrow *arrow1 = NULL;
    TArrow *arrow2 = NULL;

    TH1F* h_frame_2X1[2] = {NULL,NULL};
    TH1F* h_frame_3X1[3] = {NULL};
    TH1F* h_frame_2X4[8] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TH1F* h_frame_3X3[9] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TCanvas* c_2X1 = NULL;
    TCanvas* c_3X1 = NULL;
    TCanvas* c_2X4 = NULL;
    TCanvas* c_3X3 = NULL;
    TCanvas* c_correlate = NULL;
    TCanvas* c_single_pid = NULL;
    TCanvas* c_1X1_v2 = NULL;
    TCanvas* c_1X1_dNdpt = NULL;
    TH1F* h_frame_single_pid = NULL;
    TGraph* tg_label_plot[3][N_masses];

    ROOT::Fit::FitResult resultC;

    Tblastwave_yield_and_v2 bw_ana;
    TFile* outputfile;
    vector< vector<TGTransientFrame*> > TransientFrame_Set;
    vector< vector<TGVerticalFrame*> >  fVSetClicked;
    vector< vector<TGComboBox*> >       ComboEnergy_PID;
    vector< vector<TGComboBox*> >       ComboCentrality_PID;
    vector< vector<TGTextButton*> >     Set_energy;
    vector< vector<TGTextButton*> >     Set_centrality;
    vector <TString> vec_tgae_id_v2;
    vector <TString> vec_tgae_id_dNdpt;
    vector <TString> vec_tgae_id_v2_fit;
    vector <TString> vec_tgae_id_dNdpt_fit;

    Int_t i_transient_frame[2][N_masses_all] = {0,0};
    Double_t fit_params[4];
    TFile* RootFileFitParams;
    TString cent_upper, cent_lower;
    TString combo_energy;

public:
    TBlastWaveGUI();
    virtual ~TBlastWaveGUI();
    void CloseWindow();
    void DoText(const char *text);
    void DoSlider();
    void DoMinimize_ana();
    void DoMinimize();
    void StopMinimize();
    void TakeParamsFromMC();
    void TakeParamsFromSet();
    void Plot_curves_ana(Double_t T_BW,Double_t  Rho0_BW,Double_t  Rho2_BW,Double_t  RxOverRy_BW);
    void MakePlotv2();
    void MakePlotdNdpT();
    void CalcMaxPtLimits();
    void DrawEllipse();
    void DoSave();
    void DoSetClicked(Int_t i_particle, Int_t i_type);
    void DoSetTextButton();
    Int_t DoCentralityCombo(Int_t i_particle, Int_t i_type);
    void SetIndex(Int_t i_particle, Int_t i_type);
    void PlotData();
    Int_t DoSetSingleParticle(Int_t i_particle, Int_t i_type);
    void DoFillTgaeID();
    void WriteParams();
    ClassDef(TBlastWaveGUI, 0)
};

//______________________________________________________________________________
TBlastWaveGUI::TBlastWaveGUI() : TGMainFrame(gClient->GetRoot(), 100, 100)
{
    //--------------------------------------------------------------------
    cout << "v2_slider started" << endl;

    TGaxis::SetMaxDigits(3);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    var_test = 2.6;

    outputfile = new TFile("./out_BW.root","RECREATE");
    //--------------------------------------------------------------------


    //--------------------------------------------------------------------
    min_max_pT_range_pid[0][0] = 0.4; // pi
    min_max_pT_range_pid[1][0] = 1.0;
    min_max_pT_range_pid[0][1] = 0.4; // pi
    min_max_pT_range_pid[1][1] = 1.0;
    min_max_pT_range_pid[0][2] = 0.15; // K
    min_max_pT_range_pid[1][2]  = 1.5;
    min_max_pT_range_pid[0][3] = 0.15; // K
    min_max_pT_range_pid[1][3] = 1.5;
    min_max_pT_range_pid[0][4] = 0.22; // p
    min_max_pT_range_pid[1][4] = 1.8;
    min_max_pT_range_pid[0][5] = 0.22; // p
    min_max_pT_range_pid[1][5] = 1.8;
    min_max_pT_range_pid[0][6] = 0.5; // phi
    min_max_pT_range_pid[1][6] = 1.8;
    min_max_pT_range_pid[0][7] = 0.3; // Xi
    min_max_pT_range_pid[1][7] = 2.0;
    min_max_pT_range_pid[0][8] = 0.3; // Xi
    min_max_pT_range_pid[1][8] = 2.0;
    min_max_pT_range_pid[0][9] = 1.1; // Omega
    min_max_pT_range_pid[1][9] = 2.6;
    min_max_pT_range_pid[0][10] = 1.1; // Omega
    min_max_pT_range_pid[1][10] = 2.6;
    min_max_pT_range_pid[0][11] = 0.3; // Lambda
    min_max_pT_range_pid[1][11] = 2.0;
    min_max_pT_range_pid[0][12] = 0.3; // Lambda
    min_max_pT_range_pid[1][12] = 2.0;
    min_max_pT_range_pid[0][13] = 0.3; // K0S
    min_max_pT_range_pid[1][13] = 1.5;
    min_max_pT_range_pid[0][14] = 0.92; // D0
    min_max_pT_range_pid[1][14] = 3.6;
    min_max_pT_range_pid[0][15] = 0.22; // J/Psi
    min_max_pT_range_pid[1][15] = 4.0;
    min_max_pT_range_pid[0][16] = 0.1; // Upsilon
    min_max_pT_range_pid[1][16] = 10.0;
    min_max_pT_range_pid[0][17] = 0.5; // d
    min_max_pT_range_pid[1][17] = 2.6;
    min_max_pT_range_pid[0][18] = 0.5; // d
    min_max_pT_range_pid[1][18] = 2.6;
    min_max_pT_range_pid[0][19] = 0.3; // He3
    min_max_pT_range_pid[1][19] = 3.6;
    min_max_pT_range_pid[0][20] = 0.3; // He3
    min_max_pT_range_pid[1][20] = 3.6;
    min_max_pT_range_pid[0][21] = 0.3; // t
    min_max_pT_range_pid[1][21] = 3.6;
    //--------------------------------------------------------------------

    tg_leg = new TGraph();

    for(int i_mass = 0; i_mass < N_masses_all; ++i_mass)
    {
        tg_v2_BW_ana_pid[i_mass]    = NULL;
        tg_dNdpT_BW_ana_pid[i_mass] = NULL;

        tg_v2_BW_ana_pid_min[i_mass]    = NULL;
        tg_dNdpT_BW_ana_pid_min[i_mass] = NULL;

        tg_v2_BW_ana_pid_plot[i_mass]    = NULL;
        tg_v2_BW_ana_pid_plot_range[i_mass]    = NULL;
        tg_dNdpT_BW_ana_pid_plot[i_mass] = NULL;
        tg_dNdpT_BW_ana_pid_plot_range[i_mass] = NULL;

        TLatex_legend_dNdpT[i_mass] = NULL;
        TLatex_legend_dNdpT_plot[i_mass] = NULL;

        tg_label_plot[0][i_mass] = new TGraph();
        tg_label_plot[1][i_mass] = new TGraph();
        tg_label_plot[2][i_mass] = new TGraph();
    }


    //------------------------------------------------------------
    vec_Hframe.resize(5);
    vec_slider.resize(5);
    vec_LayoutHints.resize(5);
    vec_TextEntry.resize(5);
    vec_TextBuffer.resize(5);
    vec_TGLabel.resize(5);

    h_dummy = new TH1D("h_dummy","h_dummy",100,0,15);
    h_dummy ->GetXaxis()->SetRangeUser(0,15);
    h_dummy ->GetYaxis()->SetRangeUser(-0.2,1.1);
    h_dummy ->GetXaxis()->CenterTitle();
    h_dummy ->GetYaxis()->CenterTitle();
    h_dummy ->SetTitle("");
    h_dummy ->GetXaxis()->SetTitleOffset(1.2);
    h_dummy ->GetYaxis()->SetTitleOffset(1.1);
    h_dummy ->GetXaxis()->SetLabelSize(0.06);
    h_dummy ->GetYaxis()->SetLabelSize(0.06);
    h_dummy ->GetXaxis()->SetTitleSize(0.06);
    h_dummy ->GetYaxis()->SetTitleSize(0.06);
    h_dummy ->GetXaxis()->SetNdivisions(505,'N');
    h_dummy ->GetYaxis()->SetNdivisions(505,'N');
    h_dummy ->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_dummy ->GetYaxis()->SetTitle("v_{2}");

    h_dummy_dNdpT = new TH1D("h_dummy_dNdpT","h_dummy_dNdpT",100,0,15);
    h_dummy_dNdpT ->GetXaxis()->SetRangeUser(0,6.2);
    h_dummy_dNdpT ->GetYaxis()->SetRangeUser(-0.2,2.2);
    h_dummy_dNdpT ->GetXaxis()->CenterTitle();
    h_dummy_dNdpT ->GetYaxis()->CenterTitle();
    h_dummy_dNdpT ->SetTitle("");
    h_dummy_dNdpT ->GetXaxis()->SetTitleOffset(1.2);
    h_dummy_dNdpT ->GetYaxis()->SetTitleOffset(1.0);
    h_dummy_dNdpT ->GetXaxis()->SetLabelSize(0.06);
    h_dummy_dNdpT ->GetYaxis()->SetLabelSize(0.06);
    h_dummy_dNdpT ->GetXaxis()->SetTitleSize(0.06);
    h_dummy_dNdpT ->GetYaxis()->SetTitleSize(0.06);
    h_dummy_dNdpT ->GetXaxis()->SetNdivisions(505,'N');
    h_dummy_dNdpT ->GetYaxis()->SetNdivisions(505,'N');
    h_dummy_dNdpT ->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_dummy_dNdpT ->GetYaxis()->SetTitle("1/p_{T} dN/dp_{T} (GeV/c)^{-2}");

    TString vec_label[5] = {"T (GeV)","rho0","rhoa","Rx","fboost"};


    init_data();
    init_pT_spectra_data();
    load_data();

    //make_5_60_spectra();
    Init_v2_Mathematica();

    /*
    inputfiles[0] = TFile::Open("./Data/merge_out_v2_boost.root");
    inputfiles[1] = TFile::Open("./Data/Merge_v2_dNdpT_pT_d.root"); // deuterons

    for(Int_t i_file = 0; i_file < 2; i_file++)
    {
        for(Int_t i_R_x = 0; i_R_x < 9; i_R_x++)
        {
            for(Int_t i_fboost = 0; i_fboost < 9; i_fboost++)
            {
                printf("i_R_x: %d, i_fboost: %d \n",i_R_x,i_fboost);
                arr_lists[i_file][i_R_x][i_fboost] = (TList*)inputfiles[i_file]->Get(Form("list_BW_Rx%d_fb%d",i_R_x,i_fboost));
            }
        }
    }


    // Determine integration range of dNdpT for data
    printf("Determine integration range of dNdpT for data \n");
    for(Int_t i_mass = 0; i_mass < N_masses; i_mass++)
    {
        Double_t x_val_data_first, y_val_data_first, x_val_data_last, y_val_data_last, x_err_low_data, x_err_high_data;
        tgae_dN_dpT_mesons_data[i_mass] ->GetPoint(0,x_val_data_first,y_val_data_first);
        x_err_low_data = tgae_dN_dpT_mesons_data[i_mass] ->GetErrorXlow(0);
        tgae_dN_dpT_mesons_data[i_mass] ->GetPoint(tgae_dN_dpT_mesons_data[i_mass] ->GetN()-1,x_val_data_last,y_val_data_last);
        x_err_high_data = tgae_dN_dpT_mesons_data[i_mass] ->GetErrorXhigh(tgae_dN_dpT_mesons_data[i_mass] ->GetN()-1);

        integration_range_pid[i_mass][0] = x_val_data_first - x_err_low_data;
        integration_range_pid[i_mass][1] = x_val_data_last  + x_err_high_data;
    }


    Int_t N_masses_A[2] = {8,1};
    for(Int_t i_file = 0; i_file < 2; i_file++)
    {
        for(Int_t i_R_x = 0; i_R_x < 9; i_R_x++)
        {
            printf("i_R_x: %d \n",i_R_x);
            for(Int_t i_fboost = 0; i_fboost < 9; i_fboost++)
            {
                for(Int_t i_Temp = 0; i_Temp < 9; i_Temp++)
                {
                    //printf("i_Temp: %d \n",i_Temp);
                    for(Int_t i_rho_0 = 0; i_rho_0 < 9; i_rho_0++)
                    {
                        //printf("i_rho_0: %d \n",i_rho_0);
                        for(Int_t i_rho_a = 0; i_rho_a < 9; i_rho_a++)
                        {
                            for(Int_t i_mass = 0; i_mass < N_masses_A[i_file]; i_mass++)
                            {
                                Int_t i_mass_use = i_mass;
                                if(i_file > 0) i_mass_use = i_mass + N_masses_A[i_file-1];
                                Int_t i_N = i_mass + N_masses_A[i_file]*i_rho_a + N_masses_A[i_file]*9*i_rho_0 + N_masses_A[i_file]*9*9*i_Temp;
                                Int_t i_pos_v2    = i_N*2;
                                Int_t i_pos_dNdpT = i_N*2 + 1;
                                tp_v2_vs_pT_mesons[i_mass_use][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] = (TProfile*)arr_lists[i_file][i_R_x][i_fboost] ->At(i_pos_v2);
                                h_dN_dpT_mesons[i_mass_use][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]    = (TH1F*)arr_lists[i_file][i_R_x][i_fboost]     ->At(i_pos_dNdpT);

                                // integrate only in range of available data
                                Int_t start_integrate_bin = h_dN_dpT_mesons[i_mass_use][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->FindBin(integration_range_pid[i_mass_use][0]);
                                Int_t stop_integrate_bin  = h_dN_dpT_mesons[i_mass_use][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->FindBin(integration_range_pid[i_mass_use][1]);
                                Double_t integral = h_dN_dpT_mesons[i_mass_use][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->Integral(start_integrate_bin,stop_integrate_bin,"width");
                                //Double_t integral = h_dN_dpT_mesons[i_mass_use][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->Integral(1,-1,"width");
                                if(integral <= 0.0) continue;
                                h_dN_dpT_mesons[i_mass_use][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->Scale(1.0/integral);
                            }
                        }
                    }
                }
            }
            }
    }
    printf("All spectra normalized \n");
    */

    char buf[32];
    SetCleanup(kDeepCleanup);

    // output
    h_fit_params = new TH1D("fit_params", "fit_params", 4, 1,10);
    h_min_val_pT = new TH1D("min_val_pT", "min_val_pT", 22,1,10);
    h_max_val_pT = new TH1D("max_val_pT", "max_val_pT", 22,1,10);
    //------------------------------------------------------------



    //------------------------------------------------------------
    // Create an embedded canvas and add to the main frame, centered in x and y
    // and with 30 pixel margins all around
    fCanvas = new TRootEmbeddedCanvas("Canvas", this, 800, 700);
    fLcan = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 10);
    AddFrame(fCanvas, fLcan);
    fCanvas->GetCanvas()->SetFillColor(10);
    fCanvas->GetCanvas()->SetFrameFillColor(10);
    fCanvas->GetCanvas()->SetBorderMode(0);
    fCanvas->GetCanvas()->SetGrid();
    fCanvas->GetCanvas()->SetLogy(0);
    fCanvas->GetCanvas()->SetTopMargin(0.05);
    fCanvas->GetCanvas()->SetBottomMargin(0.17);
    fCanvas->GetCanvas()->SetRightMargin(0.02);
    fCanvas->GetCanvas()->SetLeftMargin(0.13);
    fCanvas->GetCanvas()->SetTicks(1,1);
    fCanvas->GetCanvas()->SetGrid(0,0);

    SetWindowName("Elliptic flow data");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();
    this->Resize(1200,800);

    //----------------------------------------------------
    // Slider start XA
    //Int_t start_pos_slider[5] = {5,7,3,5,0};
    //Int_t start_pos_slider[5] = {5,4,3,5,0};
    //Int_t start_pos_slider[5] = {4,4,2,6,0}; // overall OK
    //Int_t start_pos_slider[5] = {4,7,3,6,0}; // very good for v2 only
    Int_t start_pos_slider[5] = {3,5,3,7,0}; // overall quite good

    // Add a frame for the sliders
    FrameA = new TGMainFrame(gClient->GetRoot(), 400, 100);
    // Add a slider
    for(Int_t i_param = 0; i_param < 5; i_param++)
    {
        GroupSlider[i_param] = new TGGroupFrame(FrameA,"Set " + vec_label[i_param],kHorizontalFrame);
        vec_Hframe[i_param] = new TGHorizontalFrame(GroupSlider[i_param], 0, 0);
        vec_slider[i_param] = new TGHSlider(GroupSlider[i_param],250,kSlider1|kScaleDownRight,1);
        vec_slider[i_param]->Connect("PositionChanged(Int_t)", "TBlastWaveGUI",this, "DoSlider()");
        vec_slider[i_param]->SetRange(0,8);
        vec_slider[i_param]->SetPosition(start_pos_slider[i_param]);
        vec_LayoutHints[i_param] = new TGLayoutHints(kLHintsBottom | kLHintsCenterY, 5, 1, 1, 1); // handles size of slider and number box
        GroupSlider[i_param]->AddFrame(vec_slider[i_param], vec_LayoutHints[i_param]);
        FrameA->AddFrame(vec_Hframe[i_param], vec_LayoutHints[i_param]);


        // Add number text box
        vec_TextEntry[i_param] = new TGTextEntry(GroupSlider[i_param], vec_TextBuffer[i_param] = new TGTextBuffer(5), HId1);
        //vec_TextEntry[i_param] ->SetDefaultSize(0,0);
        vec_TextEntry[i_param] ->SetTextColor(kRed);
        vec_TextEntry[i_param]->SetToolTipText("Minimum (left) Value of Slider");
        vec_TextBuffer[i_param]->AddText(0, "0.0");
        vec_TextEntry[i_param]->Connect("TextChanged(char*)", "TBlastWaveGUI", this,"DoText(char*)");
        //vec_Hframe[i_param]->Resize(100, 25);
        GroupSlider[i_param]->AddFrame(vec_TextEntry[i_param], vec_LayoutHints[i_param]);
        FrameA->AddFrame(GroupSlider[i_param], vec_LayoutHints[i_param]);
        //vec_TGLabel[i_param] = new TGLabel(vec_Hframe[i_param], vec_label[i_param].Data());
        //Group[i_param]->AddFrame(vec_TGLabel[i_param], new TGLayoutHints(kLHintsBottom|kLHintsCenterY,1,1,1,1));

        vec_Hframe[i_param] ->SetHeight(1);
        //vec_Hframe[i_param] ->DrawBorder();
    }
    //----------------------------------------------------
    FrameA->SetWindowName("Sliders for Elliptic flow data");
    FrameA->MapSubwindows();
    FrameA->Resize(GetDefaultSize());
    FrameA->MapWindow();
    FrameA->Resize(435,300);
    //------------------------------------------------------------



    //------------------------------------------------------------
    FrameB = new TGMainFrame(gClient->GetRoot(), 400, 100);
    FrameB ->SetWindowName("Transverse momentum spectra");
    FrameB ->MapSubwindows();
    FrameB ->Resize(GetDefaultSize());
    FrameB ->MapWindow();
    FrameB->Resize(1300,700);


    // Create an embedded canvas and add to the main frame, centered in x and y
    // and with 30 pixel margins all around
    fCanvasB = new TRootEmbeddedCanvas("CanvasB", FrameB, 1200, 700);
    fLcanC = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 10);
    FrameB ->AddFrame(fCanvasB, fLcanC);
    fCanvasB->GetCanvas()->Divide(3,3);
    for(Int_t iPad = 1; iPad <= N_masses; iPad++)
    {
         fCanvasB->GetCanvas()->cd(iPad)->SetLogy(1);
         fCanvasB->GetCanvas()->cd(iPad);
         fCanvasB->GetCanvas()->cd(iPad)->SetFillColor(10);
         fCanvasB->GetCanvas()->cd(iPad)->SetFrameFillColor(10);
         fCanvasB->GetCanvas()->cd(iPad)->SetBorderMode(0);
         fCanvasB->GetCanvas()->cd(iPad)->SetGrid();
         fCanvasB->GetCanvas()->cd(iPad)->SetLogy(0);
         fCanvasB->GetCanvas()->cd(iPad)->SetTopMargin(0.05);
         fCanvasB->GetCanvas()->cd(iPad)->SetBottomMargin(0.17);
         fCanvasB->GetCanvas()->cd(iPad)->SetRightMargin(0.02);
         fCanvasB->GetCanvas()->cd(iPad)->SetLeftMargin(0.165);
         fCanvasB->GetCanvas()->cd(iPad)->SetTicks(1,1);
         fCanvasB->GetCanvas()->cd(iPad)->SetGrid(0,0);
    }

    // Set a name to the main frame
    FrameB->SetWindowName("Transverse momentum spectra");

    // Map all subwindows of main frame
    FrameB->MapSubwindows();

    // Map main frame
    FrameB->MapWindow();
    FrameB ->Move(400,10);

    for(Int_t iPad = 1; iPad <= N_masses; iPad++)
    {
        fCanvasB ->GetCanvas()->cd(iPad);
        h_dummy  ->Draw();
        fCanvasB ->GetCanvas()->cd(iPad)->Modified();
        fCanvasB ->GetCanvas()->cd(iPad)->Update();
    }
    //------------------------------------------------------------

    TGGC myGC = *gClient->GetResourcePool()->GetFrameGC();
    TGFont *myfont = gClient->GetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");

    //------------------------------------------------------------
    // Create horizontal splitter
    FrameD = new TGMainFrame(gClient->GetRoot(), 400, 100);
    FrameD ->SetWindowName("Buttons");

    //--------------
    // A horizontal frame
    //hframeD1  = new TGHorizontalFrame(FrameD,200,100);
    fGroupFrames[0] = new TGGroupFrame(FrameD, new TGString("Basics"),kHorizontalFrame|kRaisedFrame);

    // exit button
    Button_exit = new TGTextButton(fGroupFrames[0], "&Exit ","gApplication->Terminate(0)");
    fGroupFrames[0]->AddFrame(Button_exit, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // save button
    Button_save = new TGTextButton(fGroupFrames[0], "&Save ",10);
    Button_save->Connect("Clicked()", "TBlastWaveGUI", this, "DoSave()");
    fGroupFrames[0]->AddFrame(Button_save, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    FrameD ->AddFrame(fGroupFrames[0], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------
    //hframeD5  = new TGHorizontalFrame(FrameD,200,100);

    Button_take_params_MC_to_ana = new TGTextButton(fGroupFrames[0], "Use MC params for Ana",10);
    Button_take_params_MC_to_ana->Connect("Clicked()", "TBlastWaveGUI", this, "TakeParamsFromMC()");
    fGroupFrames[0]->AddFrame(Button_take_params_MC_to_ana, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    Button_take_params_Set_to_ana = new TGTextButton(fGroupFrames[0], "Use SET params for Ana",10);
    Button_take_params_Set_to_ana->Connect("Clicked()", "TBlastWaveGUI", this, "TakeParamsFromSet()");
    fGroupFrames[0]->AddFrame(Button_take_params_Set_to_ana, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // FrameD ->AddFrame(hframeD5, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

    //--------------

    // hframeD6  = new TGHorizontalFrame(FrameD,200,100);

    Button_make_plot_v2 = new TGTextButton(fGroupFrames[0], "Make plot v2",10);
    Button_make_plot_v2->Connect("Clicked()", "TBlastWaveGUI", this, "MakePlotv2()");
    fGroupFrames[0]->AddFrame(Button_make_plot_v2, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    Button_make_plot_dNdpT = new TGTextButton(fGroupFrames[0], "Make plot dN/dpT",10);
    Button_make_plot_dNdpT->Connect("Clicked()", "TBlastWaveGUI", this, "MakePlotdNdpT()");
    fGroupFrames[0]->AddFrame(Button_make_plot_dNdpT, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // FrameD ->AddFrame(hframeD6, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

    //--------------

    // new plot data button

    Button_plot_data= new TGTextButton(fGroupFrames[0], "Plot data",10);
    Button_plot_data->Connect("Clicked()", "TBlastWaveGUI", this, "PlotData()");
    fGroupFrames[0]->AddFrame(Button_plot_data, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    h_dummy_plot_data = new TH1D("h_dummy_plot_data","h_dummy_plot_data" , 200,0,20);
    //--------------

    // write params button
    Add_output_file_name = new TGTextEntry(fGroupFrames[0], "OutputFileName");
    fGroupFrames[0]->AddFrame(Add_output_file_name, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    Button_write_params = new TGTextButton(fGroupFrames[0], "Write params",10);
    Button_write_params ->Connect("Clicked()", "TBlastWaveGUI", this, "WriteParams()");
    fGroupFrames[0]->AddFrame(Button_write_params, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    /*
    //--------------
    // A horizontal frame
    hVframeD3  = new TGVerticalFrame(FrameD,200,100);
    fHint2 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 5, 5,  5, 10);
    fHProg1 = new TGHProgressBar(hVframeD3, 400);
    fHProg1->ShowPosition();
    hVframeD3->AddFrame(fHProg1, fHint2);
    FrameD ->AddFrame(hVframeD3, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    fHProg1->Reset();
    //--------------



    //--------------
    // A horizontal frame

    printf("Add minimize buttons \n");
    hframeD2  = new TGHorizontalFrame(FrameD,200,100);

    // exit button
    Button_minimize = new TGTextButton(hframeD2, "Minimize MC",10);
    Button_minimize->Connect("Clicked()", "TBlastWaveGUI", this, "DoMinimize()");

    Pixel_t red;
    gClient->GetColorByName("red", red);
    Button_minimize->ChangeBackground(red);

    hframeD2->AddFrame(Button_minimize, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // save button
    Button_stop_minimize = new TGTextButton(hframeD2, "Stop minimize MC",10);
    Button_stop_minimize->Connect("Clicked()", "TBlastWaveGUI", this, "StopMinimize()");
    hframeD2->AddFrame(Button_stop_minimize, new TGLayoutHints(kLHintsCenterX,5,5,3,4));


    FrameD ->AddFrame(hframeD2, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------


    */
    //--------------
    // A horizontal frame
    printf("Add progress bar \n");
    hVframeD4  = new TGVerticalFrame(FrameD,200,100);
    fHProg2 = new TGHProgressBar(hVframeD4, 400);

    fHProg2->SetBarColor("lightblue");
    fHProg2->SetRange(0.0,1000.0);
    fHProg2->ShowPosition(kTRUE, kFALSE, "%.0f calls");

    //fHProg2->ShowPosition();
    hVframeD4->AddFrame(fHProg2, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 5, 5,  5, 10));
    FrameD ->AddFrame(hVframeD4, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    fHProg2->Reset();


    //--------------


    //--------------
    // A horizontal frame
    hframeD2b  = new TGHorizontalFrame(FrameD,200,100);
    Pixel_t red;
    gClient->GetColorByName("red", red);
    // exit button
    Button_minimize_ana = new TGTextButton(hframeD2b, "Minimize Ana",10);
    Button_minimize_ana->Connect("Clicked()", "TBlastWaveGUI", this, "DoMinimize_ana()");

    Button_minimize_ana->ChangeBackground(red);

    hframeD2b->AddFrame(Button_minimize_ana, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // save button
    Button_stop_minimize_ana = new TGTextButton(hframeD2b, "Stop minimize Ana",10);
    Button_stop_minimize_ana->Connect("Clicked()", "TBlastWaveGUI", this, "StopMinimize()");
    hframeD2b->AddFrame(Button_stop_minimize_ana, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    Button_draw_ellipse_ana = new TGTextButton(hframeD2b, "Draw ellipse",10);
    Button_draw_ellipse_ana ->Connect("Clicked()", "TBlastWaveGUI", this, "DrawEllipse()");
    hframeD2b->AddFrame(Button_draw_ellipse_ana, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    fCombo = new TGComboBox(hframeD2b, 88);
    hframeD2b->AddFrame(fCombo, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    fCombo->AddEntry("fos1", 1); // used for fits as in paper -> freeze-out-hypersurface as described in paper
    fCombo->AddEntry("don't use", 2); // test freeze-out-hypersurface, not physical
    fCombo->AddEntry("boost", 3); // same as in Monte Carlo, no hypersurface, gives identical results to Alex MC
    fCombo->Resize(100, 20);
    fCombo->Select(1);


    FrameD ->AddFrame(hframeD2b, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    //--------------
    ULong_t green;
    gClient->GetColorByName("green", green);
    printf("Add particle check boxes \n");
    fGroupFrames[1] = new TGGroupFrame(FrameD, new TGString("PID fit"),kVerticalFrame|kRaisedFrame);
    fGroupFrame_PID_fit[0] = new TGGroupFrame(fGroupFrames[1],new TGString("v2"),  kVerticalFrame);
    fGroupFrame_PID_fit[1] = new TGGroupFrame(fGroupFrames[1],new TGString("dNdpT"),  kVerticalFrame);
    fGroupFrames[1]->AddFrame(fGroupFrame_PID_fit[0]);
    fGroupFrames[1]->AddFrame(fGroupFrame_PID_fit[1]);
    fCompositeFrame_pid[0] = new TGCompositeFrame(fGroupFrame_PID_fit[0],60,20, kHorizontalFrame);
    fCompositeFrame_pid[1] = new TGCompositeFrame(fGroupFrame_PID_fit[0],60,20, kHorizontalFrame);
    fCompositeFrame_pid[2] = new TGCompositeFrame(fGroupFrame_PID_fit[1],60,20, kHorizontalFrame);
    fCompositeFrame_pid[3] = new TGCompositeFrame(fGroupFrame_PID_fit[1],60,20, kHorizontalFrame);
    fGroupFrame_PID_fit[0]->AddFrame(fCompositeFrame_pid[0], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    fGroupFrame_PID_fit[0]->AddFrame(fCompositeFrame_pid[1], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    fGroupFrame_PID_fit[1]->AddFrame(fCompositeFrame_pid[2], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    fGroupFrame_PID_fit[1]->AddFrame(fCompositeFrame_pid[3], new TGLayoutHints(kLHintsCenterX,2,2,2,2));

    fLCheckBox = new TGLayoutHints(kLHintsTop | kLHintsLeft,0, 0, 5, 0);
    fLGroupFrames_particles = new TGLayoutHints(kLHintsCenterX,2,2,2,2);
    for(Int_t i_particle = 0; i_particle < 12; i_particle++)    // 20= N_masses
    {
        fGroupFrames_particles[0][i_particle] = new TGGroupFrame(fCompositeFrame_pid[0], new TGString(label_full_pid_spectra[i_particle].Data()),  kVerticalFrame);
        fGroupFrames_particles[2][i_particle] = new TGGroupFrame(fCompositeFrame_pid[2], new TGString(label_full_pid_spectra[i_particle].Data()),  kVerticalFrame);
        //fGroupFrames_particles[0][i_particle]->ChangeBackground(green);

        fCompositeFrame_pid[0]->AddFrame(fGroupFrames_particles[0][i_particle],  fLGroupFrames_particles);
        fCompositeFrame_pid[2]->AddFrame(fGroupFrames_particles[2][i_particle],  fLGroupFrames_particles);

        fCheckBox_pid[i_particle]  = new TGCheckButton(fGroupFrames_particles[0][i_particle], new TGHotString(label_full_pid_spectra[i_particle].Data()), -1);
        fCheckBox_pid_fit_dNdpt[i_particle]  = new TGCheckButton(fGroupFrames_particles[2][i_particle], new TGHotString(label_full_pid_spectra[i_particle].Data()), -1);

        fCheckBox_pid_set[0][i_particle] = new TGCheckButton(fGroupFrames_particles[0][i_particle],"set" );
        fCheckBox_pid_set[2][i_particle] = new TGCheckButton(fGroupFrames_particles[2][i_particle],"set" );

        fCheckBox_pid_set[0][i_particle]->Connect("Clicked()", "TBlastWaveGUI", this, Form("DoSetClicked(Int_t =%d, %d)", i_particle,0));
        fCheckBox_pid_set[2][i_particle]->Connect("Clicked()", "TBlastWaveGUI", this, Form("DoSetClicked(Int_t =%d, %d)", i_particle,1));

        
        //fCheckBox_pid[i_particle] ->Connect("Clicked()", "TBlastWaveGUI", this, "DoSlider()");
        fGroupFrames_particles[0][i_particle]->AddFrame(fCheckBox_pid[i_particle], fLCheckBox);
        fGroupFrames_particles[2][i_particle]->AddFrame(fCheckBox_pid_fit_dNdpt[i_particle], fLCheckBox);
        fGroupFrames_particles[0][i_particle]->AddFrame(fCheckBox_pid_set[0][i_particle], fLCheckBox);
        fGroupFrames_particles[2][i_particle]->AddFrame(fCheckBox_pid_set[2][i_particle], fLCheckBox);
    }
    for(Int_t i_particle = 12; i_particle < N_masses_all; i_particle++)    // 20= N_masses
    {
        fGroupFrames_particles[1][i_particle] = new TGGroupFrame(fCompositeFrame_pid[1], new TGString(label_full_pid_spectra[i_particle].Data()),  kVerticalFrame);
        fGroupFrames_particles[3][i_particle] = new TGGroupFrame(fCompositeFrame_pid[3], new TGString(label_full_pid_spectra[i_particle].Data()),  kVerticalFrame);

        fCompositeFrame_pid[1]->AddFrame(fGroupFrames_particles[1][i_particle],  fLGroupFrames_particles);
        fCompositeFrame_pid[3]->AddFrame(fGroupFrames_particles[3][i_particle],  fLGroupFrames_particles);

        fCheckBox_pid[i_particle]            = new TGCheckButton(fGroupFrames_particles[1][i_particle], new TGHotString(label_full_pid_spectra[i_particle].Data()), -1);
        fCheckBox_pid_fit_dNdpt[i_particle]  = new TGCheckButton(fGroupFrames_particles[3][i_particle], new TGHotString(label_full_pid_spectra[i_particle].Data()), -1);

        fCheckBox_pid_set[1][i_particle] = new TGCheckButton(fGroupFrames_particles[1][i_particle],"set" );
        fCheckBox_pid_set[3][i_particle] = new TGCheckButton(fGroupFrames_particles[3][i_particle],"set" );

        fCheckBox_pid_set[1][i_particle]->Connect("Clicked()", "TBlastWaveGUI", this, Form("DoSetClicked(Int_t= %d, %d)", i_particle,0));
        fCheckBox_pid_set[3][i_particle]->Connect("Clicked()", "TBlastWaveGUI", this, Form("DoSetClicked(Int_t= %d, %d)", i_particle,1));

        //fCheckBox_pid[i_particle] ->Connect("Clicked()", "TBlastWaveGUI", this, "DoSlider()");
        fGroupFrames_particles[1][i_particle]->AddFrame(fCheckBox_pid[i_particle], fLCheckBox);
        fGroupFrames_particles[3][i_particle]->AddFrame(fCheckBox_pid_fit_dNdpt[i_particle], fLCheckBox);
        fGroupFrames_particles[1][i_particle]->AddFrame(fCheckBox_pid_set[1][i_particle], fLCheckBox);
        fGroupFrames_particles[3][i_particle]->AddFrame(fCheckBox_pid_set[3][i_particle], fLCheckBox);
    }
    FrameD ->AddFrame(fGroupFrames[1], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------


    //--------------
    TransientFrame_Set.resize(2); // v2, dNdpT
    fVSetClicked.resize(2);
    Set_energy.resize(2);
    Set_centrality.resize(2);
    ComboCentrality_PID.resize(2);
    ComboEnergy_PID.resize(2);
    for(Int_t i_type = 0; i_type < 2; i_type++)
    {
        TransientFrame_Set[i_type].resize(N_masses_all);
        fVSetClicked[i_type].resize(N_masses_all);
        ComboCentrality_PID[i_type].resize(N_masses_all);
        ComboEnergy_PID[i_type].resize(N_masses_all);
        Set_energy[i_type].resize(N_masses_all);
        Set_centrality[i_type].resize(N_masses_all);
        for(Int_t i_particle = 0; i_particle < N_masses_all; i_particle++)
        {
            TransientFrame_Set[i_type][i_particle] = NULL;
            fVSetClicked[i_type][i_particle]       = NULL;
            ComboCentrality_PID[i_type][i_particle]= NULL;
            ComboEnergy_PID[i_type][i_particle]    = NULL;
            Set_energy[i_type][i_particle]         = NULL;
            Set_centrality[i_type][i_particle]     = NULL;
        }
    }
    //--------------


    //--------------
    printf("Add particle check boxes \n");
    fGroupFrames[5] = new TGGroupFrame(FrameD, new TGString("PID plot"),kVerticalFrame|kRaisedFrame);
    fLCheckBoxB = new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 0, 5, 0);
    fGroupFrame_PID_plot[0] = new TGGroupFrame(fGroupFrames[5],new TGString("v2"),  kVerticalFrame);
    fGroupFrame_PID_plot[1] = new TGGroupFrame(fGroupFrames[5],new TGString("dNdpT"),  kVerticalFrame);
    fGroupFrames[5]->AddFrame(fGroupFrame_PID_plot[0]);
    fGroupFrames[5]->AddFrame(fGroupFrame_PID_plot[1]);
    fCompositeFrame_pid[4] = new TGCompositeFrame(fGroupFrame_PID_plot[0],60,20, kHorizontalFrame);     //dNdpt
    fCompositeFrame_pid[5] = new TGCompositeFrame(fGroupFrame_PID_plot[0],60,20, kHorizontalFrame);     //dNdpt
    fCompositeFrame_pid[6] = new TGCompositeFrame(fGroupFrame_PID_plot[1],60,20, kHorizontalFrame);     //v2
    fCompositeFrame_pid[7] = new TGCompositeFrame(fGroupFrame_PID_plot[1],60,20, kHorizontalFrame);     //v2
    fGroupFrame_PID_plot[0]->AddFrame(fCompositeFrame_pid[4], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    fGroupFrame_PID_plot[0]->AddFrame(fCompositeFrame_pid[5], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    fGroupFrame_PID_plot[1]->AddFrame(fCompositeFrame_pid[6], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    fGroupFrame_PID_plot[1]->AddFrame(fCompositeFrame_pid[7], new TGLayoutHints(kLHintsCenterX,2,2,2,2));

    fLCheckBox = new TGLayoutHints(kLHintsTop | kLHintsLeft,2, 2, 5, 0);
    fLGroupFrames_particles = new TGLayoutHints(kLHintsCenterX,2,2,2,2);


    for(Int_t i_particle = 0; i_particle < 12; i_particle++)
    {
        fCheckBox_pid_plot[i_particle]        = new TGCheckButton(fCompositeFrame_pid[4], new TGHotString(label_full_pid_spectra[i_particle].Data()), -1);
        fCheckBox_pid_plot_dNdpt[i_particle]  = new TGCheckButton(fCompositeFrame_pid[6], new TGHotString(label_full_pid_spectra[i_particle].Data()), -1);

        fCheckBox_pid_plot[i_particle] ->Connect("Clicked()", "TBlastWaveGUI", this, "DoSlider()");

        fCompositeFrame_pid[4]->AddFrame(fCheckBox_pid_plot[i_particle], fLCheckBox);
        fCompositeFrame_pid[6]->AddFrame(fCheckBox_pid_plot_dNdpt[i_particle], fLCheckBox);
    }
    for(Int_t i_particle = 12; i_particle < N_masses_all; i_particle++)
    {
        fCheckBox_pid_plot[i_particle]            = new TGCheckButton(fCompositeFrame_pid[5], new TGHotString(label_full_pid_spectra[i_particle].Data()), -1);
        fCheckBox_pid_plot_dNdpt[i_particle]      = new TGCheckButton(fCompositeFrame_pid[7], new TGHotString(label_full_pid_spectra[i_particle].Data()), -1);

        fCheckBox_pid_plot[i_particle] ->Connect("Clicked()", "TBlastWaveGUI", this, "DoSlider()");

        fCompositeFrame_pid[5]->AddFrame(fCheckBox_pid_plot[i_particle], fLCheckBox);
        fCompositeFrame_pid[7]->AddFrame(fCheckBox_pid_plot_dNdpt[i_particle], fLCheckBox);
    }

    FrameD ->AddFrame(fGroupFrames[5], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    //--------------
    //hframeD4  = new TGHorizontalFrame(FrameD,200,100);
    fGroupFrames[2] = new TGGroupFrame(FrameD, new TGString("Fit and plot"),kHorizontalFrame|kRaisedFrame);
    for(Int_t i_v2_dNdpT = 0; i_v2_dNdpT < 2; i_v2_dNdpT++)
    {
        fCheckBox_v2_dNdpT[i_v2_dNdpT]  = new TGCheckButton(fGroupFrames[2], new TGHotString(label_v2_dNdpT[i_v2_dNdpT].Data()), -1);
        fCheckBox_v2_dNdpT[i_v2_dNdpT] ->SetState(kButtonDown);
        fGroupFrames[2]->AddFrame(fCheckBox_v2_dNdpT[i_v2_dNdpT], fLCheckBox);
    }
    for(Int_t i_cb = 0; i_cb < 3; i_cb++)
    {
        fCheckBox_sel[i_cb]  = new TGCheckButton(fGroupFrames[2], new TGHotString(label_checkbox[i_cb].Data()), -1);
        fCheckBox_sel[i_cb] ->SetState(kButtonDown);
        fCheckBox_sel[i_cb] ->Connect("Clicked()", "TBlastWaveGUI", this, "DoSlider()");
        fGroupFrames[2]->AddFrame(fCheckBox_sel[i_cb], fLCheckBox);
    }
    FrameD ->AddFrame(fGroupFrames[2], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    //--------------
    TString arr_label_params_ana[4] = {"SET T","SET rho0","SET rhoa","SET Rx"};
    //hframeD5a  = new TGHorizontalFrame(FrameD,200,100);
    fGroupFrames[3] = new TGGroupFrame(FrameD, new TGString("SET parameters"),kHorizontalFrame|kRaisedFrame);
    for(Int_t i_param = 0; i_param < 4; i_param++)
    {
        hVframeD5a[i_param] = new TGVerticalFrame(fGroupFrames[3], 200,200);
        arr_NEntry_ana_params[i_param] = new TGNumberEntry(hVframeD5a[i_param], 0.0, 12,(TGNumberFormat::EStyle) 2);
        arr_NEntry_ana_params[i_param] ->SetNumStyle( TGNumberFormat::kNESRealTwo); // https://root.cern.ch/doc/master/classTGNumberFormat.html#a8a0f81aac8ac12d0461aef554c6271ad
        hVframeD5a[i_param]->AddFrame(arr_NEntry_ana_params[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));

        TString label_entry = arr_label_params_ana[i_param];
        arr_Label_NEntry_ana_params[i_param] = new TGLabel(hVframeD5a[i_param], label_entry.Data(), myGC(), myfont->GetFontStruct());
        hVframeD5a[i_param]->AddFrame(arr_Label_NEntry_ana_params[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
        fGroupFrames[3]->AddFrame(hVframeD5a[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    }
    FrameD ->AddFrame(fGroupFrames[3], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    //--------------
    // Add feed down check box
    fCheckBoxFeedDown = new TGCheckButton(fGroupFrames[2], "Include Feed-Down");
    fCheckBoxFeedDown->SetState(kButtonDown);
    fGroupFrames[2]->AddFrame(fCheckBoxFeedDown, new TGLayoutHints(kLHintsCenterX,5,5,6,4));

    //Add energy drop down
    TString ComboEnergyLabel[10] = {"7.7    GeV","11.5  GeV","14.5  GeV","19.6  GeV","27.0  GeV","39.0  GeV","62.4  GeV","200   GeV", "2760.0 GeV", "5020.0 GeV"};
    ComboEnergy = new TGComboBox(hframeD2b,200);
    for(Int_t i = 0; i < 10; ++i)
    {
        ComboEnergy->AddEntry(ComboEnergyLabel[i], i);
    }
    hframeD2b->AddFrame(ComboEnergy, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    ComboEnergy->Resize(100,20);
    ComboEnergy->Select(0);
    //-------------
    //Add centrality drop down
    TString ComboCentralityLabel[18] = {"0-5  %", "0-10 %", "0-20 %","0-30 %","0-80 %", "5-10 %","10-20 %", "20-30 %","10-40 %","20-40 %", "30-40 %", "40-50 %", "40-60 %","40-80 %", "50-60 %", "60-70 %","60-80 %", "70-80 %"};
    ComboCentrality = new TGComboBox(hframeD2b,200);
    for(Int_t index_cent = 0; index_cent < 18; index_cent++)
    {
        ComboCentrality->AddEntry(ComboCentralityLabel[index_cent], index_cent);
    }
    hframeD2b->AddFrame(ComboCentrality, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    ComboCentrality->Resize(100,20);
    ComboCentrality->Select(0);
    //-------------

    //-------------
    Set_energy_centrality = new TGTextButton(hframeD2b, "Set");
    hframeD2b->AddFrame(Set_energy_centrality, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    Set_energy_centrality->Resize(100,20);
    Set_energy_centrality->Connect("Clicked()", "TBlastWaveGUI", this, "DoSetTextButton()");
    //-------------


    FrameD ->MapSubwindows();
    FrameD ->MapWindow();
    FrameD ->Resize(1050,900); // size of frame
    FrameD ->Move(1250,850); // position of frame



    //--------------
    printf("Add limits widget \n");
    LHintsD4a = new TGLayoutHints(kLHintsCenterX,5,5,2,0);
    LHintsD4a2 = new TGLayoutHints(kLHintsCenterX,5,5,5,0);

   // TString arr_label_pid[N_masses]     = {"Pi","K","P", "Phi","Omega","D0","J/psi","Y","d"};
    TString arr_label_min_max[2] = {"min","max"};

    arr_fL1[0] = new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 0);
    arr_fL1[1] = new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 0);

    Double_t max_val_VF[2] = {300,500};

    frame_TGTransient = new TGTransientFrame(gClient->GetRoot(), FrameD, 60, 20, kHorizontalFrame);
    FrameD ->AddFrame(frame_TGTransient,LHintsD4a);

    Hframe_pT_limits = new TGHorizontalFrame(frame_TGTransient, 60, 20);
    frame_TGTransient ->AddFrame(Hframe_pT_limits);
    Vframe_pT_limits[0] = new TGVerticalFrame(Hframe_pT_limits, 200,300);
    Hframe_pT_limits->AddFrame(Vframe_pT_limits[0],new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    Vframe_pT_limits[1] = new TGVerticalFrame(Hframe_pT_limits);
    Hframe_pT_limits->AddFrame(Vframe_pT_limits[1]);

    fGroupFrames[4] = new TGGroupFrame(Vframe_pT_limits[0], new TGString("Set pT limits"),kHorizontalFrame|kRaisedFrame);
    Vframe_pT_limits[0]->AddFrame(fGroupFrames[4], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    NEntry_set_limits = new TGNumberEntry(fGroupFrames[4], 0.68, 12,(TGNumberFormat::EStyle) 1);
    NEntry_set_limits ->Connect("ValueSet(Long_t)", "TBlastWaveGUI", this, "CalcMaxPtLimits()");
    NEntry_set_limits->SetNumStyle( TGNumberFormat::kNESRealTwo); // https://root.cern.ch/doc/master/classTGNumberFormat.html#a8a0f81aac8ac12d0461aef554c6271ad
    fGroupFrames[4]->AddFrame(NEntry_set_limits, LHintsD4a);

    for(Int_t i_pid = 0; i_pid < N_masses_all; i_pid++)
    {
        pT_Range_Group[0][i_pid] = new TGGroupFrame(Vframe_pT_limits[0],label_full_pid_spectra[i_pid],kHorizontalFrame);
        if (i_pid<11)
        {
            //pT_Range_Group[0][i_pid] = new TGGroupFrame(Vframe_pT_limits[0],label_full_pid_spectra[i_pid],kHorizontalFrame);
            Vframe_pT_limits[0]->AddFrame(pT_Range_Group[0][i_pid], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
        }
        pT_Range_Group[1][i_pid] = new TGGroupFrame(Vframe_pT_limits[1],label_full_pid_spectra[i_pid],kHorizontalFrame);
        if (i_pid>10)
        {
            //pT_Range_Group[1][i_pid] = new TGGroupFrame(Vframe_pT_limits[1],label_full_pid_spectra[i_pid],kHorizontalFrame);
            Vframe_pT_limits[1]->AddFrame(pT_Range_Group[1][i_pid], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
        }
        for(Int_t i_min_max = 0; i_min_max < 2; i_min_max++)
        {
            arr_NEntry_limits[i_min_max][i_pid]   = new TGNumberEntry(pT_Range_Group[0][i_pid], min_max_pT_range_pid[i_min_max][i_pid], 12,(TGNumberFormat::EStyle) 1);
            arr_NEntry_limits_A[i_min_max][i_pid] = new TGNumberEntry(pT_Range_Group[1][i_pid], min_max_pT_range_pid[i_min_max][i_pid], 12,(TGNumberFormat::EStyle) 1);

            arr_NEntry_limits[i_min_max][i_pid] ->Connect("ValueSet(Long_t)", "TBlastWaveGUI", this, "DoSlider()");
            //arr_NEntry_limits_A[i_min_max][i_pid] ->Connect("ValueSet(Long_t)", "TBlastWaveGUI", this, "DoSlider()");

            arr_NEntry_limits[i_min_max][i_pid]  ->SetNumStyle( TGNumberFormat::kNESRealOne); // https://root.cern.ch/doc/master/classTGNumberFormat.html#a8a0f81aac8ac12d0461aef554c6271ad
            arr_NEntry_limits_A[i_min_max][i_pid]->SetNumStyle( TGNumberFormat::kNESRealOne); // https://root.cern.ch/doc/master/classTGNumberFormat.html#a8a0f81aac8ac12d0461aef554c6271ad

            pT_Range_Group[0][i_pid]->AddFrame(arr_NEntry_limits[i_min_max][i_pid], LHintsD4a);
            pT_Range_Group[1][i_pid]->AddFrame(arr_NEntry_limits_A[i_min_max][i_pid], LHintsD4a);

            TString label_entry = arr_label_min_max[i_min_max];
            arr_Label_NEntry_limits[i_min_max][i_pid]   = new TGLabel(pT_Range_Group[0][i_pid], label_entry.Data(), myGC(), myfont->GetFontStruct());
            arr_Label_NEntry_limits_A[i_min_max][i_pid] = new TGLabel(pT_Range_Group[1][i_pid], label_entry.Data(), myGC(), myfont->GetFontStruct());

            pT_Range_Group[0][i_pid]->AddFrame(arr_Label_NEntry_limits[i_min_max][i_pid], LHintsD4a2);
            pT_Range_Group[1][i_pid]->AddFrame(arr_Label_NEntry_limits_A[i_min_max][i_pid], LHintsD4a2);

        }
    }

    frame_TGTransient->MapSubwindows();
    frame_TGTransient->Resize();
    frame_TGTransient->CenterOnParent();
    frame_TGTransient->SetWindowName("pT limits");
    frame_TGTransient->MapWindow();
    frame_TGTransient->Move(270,650); // position of frame
    //--------------


    //------------------------------------------------------------



    //------------------------------------------------------------
    DoSlider();
    //------------------------------------------------------------

}
//______________________________________________________________________________
TBlastWaveGUI::~TBlastWaveGUI()
{
    // Clean up
    Cleanup();
}
//______________________________________________________________________________
void TBlastWaveGUI::CloseWindow()
{
    // Called when window is closed via the window manager.
    delete this;
}


//______________________________________________________________________________
void TBlastWaveGUI::DoSetTextButton()
{
    // Handle text button Set, selection of energy and centrality for all particles at the same time

    cout << "DoSetTextButton()" << endl;
    Int_t i_select_energy    = ComboEnergy->GetSelected();
    Int_t i_select_centrality = ComboCentrality->GetSelected();
    TString ComboEnergyLabel[10]           = {"7.7","11.5","14.5","19.6","27","39","62.4","200", "2760", "5020"};
    TString ComboEnergyLabel_float[10]     = {"7.7","11.5","14.5","19.6","27.0", "39.0","62.4", "200.0", "2760.0", "5020.0"};
    TString ComboCentralityLabel_upper[18] = {"5","10","20","30","80","10","20", "30","40","40", "40", "50", "60","80", "60", "70","80", "80"};
    TString ComboCentralityLabel_lower[18] = {"0", "0", "0", "0", "0", "5","10", "20","10","20", "30", "40", "40","40", "50", "60","60", "70"};

    ULong_t green;
    gClient->GetColorByName("green", green);
    ULong_t red;
    gClient->GetColorByName("red", red);

    for (Int_t i_particle= 0; i_particle< N_masses_all; i_particle++)
        {
            fCheckBox_pid[i_particle]           ->ChangeBackground(red);       // Check box PID fit v2
            fCheckBox_pid_fit_dNdpt[i_particle] ->ChangeBackground(red);       // Check box PID fit dNdpt
            fCheckBox_pid_plot[i_particle]      ->ChangeBackground(red);       // Check box PID plot v2
            fCheckBox_pid_plot_dNdpt[i_particle]->ChangeBackground(red);       // Check box PID plot dNdpt
            fCheckBox_pid[i_particle]           ->SetState(kButtonUp);
            fCheckBox_pid_fit_dNdpt[i_particle] ->SetState(kButtonUp);
            fCheckBox_pid_plot[i_particle]      ->SetState(kButtonUp);
            fCheckBox_pid_plot_dNdpt[i_particle]->SetState(kButtonUp);
            for (Int_t i_found = 0; i_found < (Int_t)vec_pid_energy_v2[i_particle].size();i_found++ )
            {
                if (vec_pid_energy_v2[i_particle][i_found] == ComboEnergyLabel[i_select_energy].Data() && vec_pid_cent_upper_v2[i_particle][i_found]== ComboCentralityLabel_upper[i_select_centrality].Data() && vec_pid_cent_lower_v2[i_particle][i_found]== ComboCentralityLabel_lower[i_select_centrality].Data() )
                {
                    fCheckBox_pid[i_particle] ->ChangeBackground(green);
                    fCheckBox_pid[i_particle] ->SetState(kButtonDown);
                    fCheckBox_pid_plot[i_particle]->ChangeBackground(green);
                    fCheckBox_pid_plot[i_particle] ->SetState(kButtonDown);
                }
                if (vec_pid_energy_v2[i_particle][i_found] == ComboEnergyLabel_float[i_select_energy].Data() && vec_pid_cent_upper_v2[i_particle][i_found]== ComboCentralityLabel_upper[i_select_centrality].Data() && vec_pid_cent_lower_v2[i_particle][i_found]== ComboCentralityLabel_lower[i_select_centrality].Data() )
                {
                    fCheckBox_pid[i_particle] ->ChangeBackground(green);
                    fCheckBox_pid[i_particle] ->SetState(kButtonDown);
                    fCheckBox_pid_plot[i_particle]->ChangeBackground(green);
                    fCheckBox_pid_plot[i_particle] ->SetState(kButtonDown);
                }
                
            }
            for (Int_t i_found = 0; i_found < (Int_t)vec_pid_energy_dNdpt[i_particle].size();i_found++ )
            {
                if (vec_pid_energy_dNdpt[i_particle][i_found] == ComboEnergyLabel[i_select_energy].Data() && vec_pid_cent_upper_dNdpt[i_particle][i_found]== ComboCentralityLabel_upper[i_select_centrality].Data() && vec_pid_cent_lower_dNdpt[i_particle][i_found]== ComboCentralityLabel_lower[i_select_centrality].Data() )
                {
                    fCheckBox_pid_fit_dNdpt[i_particle] ->ChangeBackground(green);
                    fCheckBox_pid_fit_dNdpt[i_particle] ->SetState(kButtonDown);
                    fCheckBox_pid_plot_dNdpt[i_particle]->ChangeBackground(green);
                    fCheckBox_pid_plot_dNdpt[i_particle] ->SetState(kButtonDown);
                }
                if (vec_pid_energy_dNdpt[i_particle][i_found] == ComboEnergyLabel_float[i_select_energy].Data() && vec_pid_cent_upper_dNdpt[i_particle][i_found]== ComboCentralityLabel_upper[i_select_centrality].Data() && vec_pid_cent_lower_dNdpt[i_particle][i_found]== ComboCentralityLabel_lower[i_select_centrality].Data() )
                {
                    fCheckBox_pid_fit_dNdpt[i_particle] ->ChangeBackground(green);
                    fCheckBox_pid_fit_dNdpt[i_particle] ->SetState(kButtonDown);
                    fCheckBox_pid_plot_dNdpt[i_particle]->ChangeBackground(green);
                    fCheckBox_pid_plot_dNdpt[i_particle] ->SetState(kButtonDown);
                }
            }
        }


    DoFillTgaeID();

}

//______________________________________________________________________________
void TBlastWaveGUI::DoFillTgaeID()
{
    // Create identifiers of the selected  particles

    TGTextLBEntry *cent_text_entry = (TGTextLBEntry *)ComboCentrality->GetSelectedEntry();
    const char *selected_centrality = cent_text_entry->GetTitle();

    TString select_centrality = selected_centrality ;
    Ssiz_t found = select_centrality.First("-");
    TSubString sub_str = select_centrality(0,found);
    cent_lower = sub_str;

    select_centrality.Replace(0, found+1, "");
    found  = select_centrality.First(" ");
    sub_str = select_centrality(0,found);
    cent_upper = sub_str;

    TGTextLBEntry *energy_text_entry = (TGTextLBEntry *)ComboEnergy->GetSelectedEntry();
    const char *selected_Energy = energy_text_entry->GetTitle();
    TString select_Energy = selected_Energy ;
    found = select_Energy.First(" ");
    TSubString sub_str_Energy = select_Energy(0,found);
    combo_energy = sub_str_Energy;

    vec_tgae_id_v2.clear();
    vec_tgae_id_dNdpt.clear();
    vec_tgae_id_v2_fit.clear();
    vec_tgae_id_dNdpt_fit.clear();
    TString tgae_id;

    for (Int_t i_particle= 0; i_particle< N_masses_all; i_particle++)
    {
        if (fCheckBox_pid_plot[i_particle]->IsDown())
        {
            tgae_id = "v2_PID_";
            tgae_id += label_full_pid_spectra[i_particle].Data();
            tgae_id += "_E_";
            tgae_id += sub_str_Energy;
            tgae_id += "_C_";
            tgae_id += cent_lower.Data();
            tgae_id += "_";
            tgae_id += cent_upper.Data();

            cout << tgae_id << endl;
            vec_tgae_id_v2.push_back(tgae_id.Copy());  // v2 plot
            tgae_id.Clear();
        }
        if (fCheckBox_pid[i_particle]->IsDown())
        {
            tgae_id = "v2_PID_";
            tgae_id += label_full_pid_spectra[i_particle].Data();
            tgae_id += "_E_";
            tgae_id += sub_str_Energy;
            tgae_id += "_C_";
            tgae_id += cent_lower.Data();
            tgae_id += "_";
            tgae_id += cent_upper.Data();

            cout << tgae_id << endl;
            vec_tgae_id_v2_fit.push_back(tgae_id.Copy());   // v2 fit
            tgae_id.Clear();
        }

        if (fCheckBox_pid_plot_dNdpt[i_particle]->IsDown())
        {
            tgae_id = "dNdpt_PID_";
            tgae_id += label_full_pid_spectra[i_particle].Data();
            tgae_id += "_E_";
            tgae_id += sub_str_Energy;
            tgae_id += "_C_";
            tgae_id += cent_lower.Data();
            tgae_id += "_";
            tgae_id += cent_upper.Data();
            cout << tgae_id << endl;
            vec_tgae_id_dNdpt.push_back(tgae_id.Copy());   // dNdpt plot
            tgae_id.Clear();
        }
        if (fCheckBox_pid_fit_dNdpt[i_particle]->IsDown())
        {
            tgae_id = "dNdpt_PID_";
            tgae_id += label_full_pid_spectra[i_particle].Data();
            tgae_id += "_E_";
            tgae_id += sub_str_Energy;
            tgae_id += "_C_";
            tgae_id += cent_lower.Data();
            tgae_id += "_";
            tgae_id += cent_upper.Data();
            cout << tgae_id << endl;
            vec_tgae_id_dNdpt_fit.push_back(tgae_id.Copy());   // dNdpt fit
            tgae_id.Clear();
        }
    }

}

//______________________________________________________________________________
void TBlastWaveGUI::DoSetClicked(Int_t i_particle, Int_t i_type)
{
    // Called when set is clicked for single particle, transient frame opens

    cout << "SetClicked()" << endl;

    if (i_transient_frame[i_type][i_particle] == 0)
    {
        TransientFrame_Set[i_type][i_particle] = new TGTransientFrame(gClient->GetRoot(), FrameD, 60, 20, kHorizontalFrame);
        fVSetClicked[i_type][i_particle]       = new TGVerticalFrame(TransientFrame_Set[i_type][i_particle]);
        ComboEnergy_PID[i_type][i_particle]    = new TGComboBox(fVSetClicked[i_type][i_particle], "Energy");
        Set_energy[i_type][i_particle]         = new TGTextButton(fVSetClicked[i_type][i_particle], "Set energy");
        ComboCentrality_PID[i_type][i_particle]= new TGComboBox(fVSetClicked[i_type][i_particle], "Centrality");
        Set_centrality[i_type][i_particle]     = new TGTextButton(fVSetClicked[i_type][i_particle], "Set centrality");
        TransientFrame_Set[i_type][i_particle] ->AddFrame(fVSetClicked[i_type][i_particle],        new TGLayoutHints(kLHintsCenterX,5,5,5,5));
        fVSetClicked[i_type][i_particle]       ->AddFrame(ComboEnergy_PID[i_type][i_particle],     new TGLayoutHints(kLHintsCenterX,5,5,3,4));
        fVSetClicked[i_type][i_particle]       ->AddFrame(Set_energy[i_type][i_particle],          new TGLayoutHints(kLHintsCenterX,5,5,3,4));
        fVSetClicked[i_type][i_particle]       ->AddFrame(ComboCentrality_PID[i_type][i_particle], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
        fVSetClicked[i_type][i_particle]       ->AddFrame(Set_centrality[i_type][i_particle],      new TGLayoutHints(kLHintsCenterX,5,5,3,4));
        Set_centrality[i_type][i_particle]     ->Connect("Clicked()", "TBlastWaveGUI", this, Form("DoSetSingleParticle(Int_t= %d, %d)", i_particle,i_type));

        ComboCentrality_PID[i_type][i_particle] ->Resize(100,20);
        ComboCentrality_PID[i_type][i_particle] ->Select(0);
        ComboEnergy_PID[i_type][i_particle]     ->Resize(100,20);
        ComboEnergy_PID[i_type][i_particle]     ->Select(0);
        TransientFrame_Set[i_type][i_particle]  ->SetWindowName(label_full_pid_spectra[i_particle].Data());

        TString ComboEnergyLabel_PID;
        if (i_type == 0)
        {
            for(Int_t i_energy_pid = 1; i_energy_pid < vec_pid_energy_v2[i_particle].size(); i_energy_pid++)
            {
                if (i_energy_pid == 1)
                {
                    ComboEnergyLabel_PID  = vec_pid_energy_v2[i_particle][i_energy_pid];
                    ComboEnergyLabel_PID += " GeV";
                    ComboEnergy_PID[i_type][i_particle]->AddEntry(ComboEnergyLabel_PID, i_energy_pid-1);
                }
                if (i_energy_pid > 1)
                {
                    if (vec_pid_energy_v2[i_particle][i_energy_pid] != vec_pid_energy_v2[i_particle][i_energy_pid-1])
                    {
                        ComboEnergyLabel_PID  = vec_pid_energy_v2[i_particle][i_energy_pid];
                        ComboEnergyLabel_PID += " GeV";
                        ComboEnergy_PID[i_type][i_particle]->AddEntry(ComboEnergyLabel_PID, i_energy_pid-1);
                    }
                }

            }
        }

        if (i_type == 1)
        {
            for(Int_t i_energy_pid = 1; i_energy_pid < vec_pid_energy_dNdpt[i_particle].size(); i_energy_pid++)
            {
                if (i_energy_pid == 1)
                {
                    ComboEnergyLabel_PID  = vec_pid_energy_dNdpt[i_particle][1];
                    ComboEnergyLabel_PID += " GeV";
                    ComboEnergy_PID[i_type][i_particle]->AddEntry(ComboEnergyLabel_PID, 0);
                }
                if (i_energy_pid > 1)
                {
                    if (vec_pid_energy_dNdpt[i_particle][i_energy_pid] != vec_pid_energy_dNdpt[i_particle][i_energy_pid-1])
                    {
                        ComboEnergyLabel_PID  = vec_pid_energy_dNdpt[i_particle][i_energy_pid];
                        ComboEnergyLabel_PID += " GeV";
                        ComboEnergy_PID[i_type][i_particle]->AddEntry(ComboEnergyLabel_PID, i_energy_pid-1);
                    }
                }
            }
        }

        //ComboEnergy_PID[i_type][i_particle]->Connect("Selected(Int_t)", "TBlastWaseGUI", this, Form("DoCentralityCombo(Int_t= %d, %d)", i_particle, i_type));
        Set_energy[i_type][i_particle] ->Connect("Clicked()", "TBlastWaveGUI", this, Form("DoCentralityCombo(Int_t= %d, %d)", i_particle,i_type));

        TransientFrame_Set[i_type][i_particle]->MapSubwindows();
        TransientFrame_Set[i_type][i_particle]->Resize(150,150);
        TransientFrame_Set[i_type][i_particle]->CenterOnParent();
        TransientFrame_Set[i_type][i_particle]->MapWindow();
    }

    if (i_transient_frame[i_type][i_particle] == 1)
    {
        TransientFrame_Set[i_type][i_particle]->UnmapWindow();
    }

    if ( i_transient_frame[i_type][i_particle]>1)
    {
        TransientFrame_Set[i_type][i_particle]->MapWindow();
        i_transient_frame[i_type][i_particle] = i_transient_frame[i_type][i_particle] -2;
    }
    i_transient_frame[i_type][i_particle] ++;
    //TransientFrame_Set[i_type][i_particle]->Move(500,650); // position of frame


}

//______________________________________________________________________________
Int_t TBlastWaveGUI::DoCentralityCombo(Int_t i_particle, Int_t i_type)
{
    // Centrality Combo is filled with available centralities at given energy

    Int_t i_entry = ComboEnergy_PID[i_type][i_particle]->GetSelected();
    if (i_entry == -1) return 0;

    cout << "DoCentralityCombo()" << endl;
    vector<TString> vec_ComboCentralityLabel_PID;
    TString ComboCentralityLabel_PID;
    vec_ComboCentralityLabel_PID.clear();

    TGTextLBEntry *filePointer = (TGTextLBEntry *)ComboEnergy_PID[i_type][i_particle]->GetSelectedEntry();
    const char *selected_energy = filePointer->GetTitle();
    TString select_energy = selected_energy ;
    Ssiz_t found = select_energy.First(" ");
    TSubString sub_str = select_energy(0,found);

    Int_t i_number  = 0;
    Int_t i_entries = ComboCentrality_PID[i_type][i_particle]->GetNumberOfEntries();
    ComboCentrality_PID[i_type][i_particle]->RemoveAll();
    if (i_type == 0)
    {
        for (Int_t i_cent_pid = 1; i_cent_pid < vec_pid_cent_upper_v2[i_particle].size(); i_cent_pid++)
        {
            if (sub_str == vec_pid_energy_v2[i_particle][i_cent_pid])
            {
                ComboCentralityLabel_PID  = vec_pid_cent_lower_v2[i_particle][i_cent_pid];
                ComboCentralityLabel_PID += "-";
                ComboCentralityLabel_PID += vec_pid_cent_upper_v2[i_particle][i_cent_pid];
                ComboCentralityLabel_PID += " %";
                if ( i_number == 0 ) ComboCentrality_PID[i_type][i_particle]->AddEntry(ComboCentralityLabel_PID, i_cent_pid-1);
                vec_ComboCentralityLabel_PID.push_back(ComboCentralityLabel_PID.Copy());
                if (i_number>0)
                {
                    if (vec_ComboCentralityLabel_PID[i_number -1] != vec_ComboCentralityLabel_PID[i_number]) ComboCentrality_PID[i_type][i_particle]->AddEntry(ComboCentralityLabel_PID, i_cent_pid-1);
                }
                ComboCentralityLabel_PID.Clear();
                i_number ++;
            }
        }

    }

    if (i_type == 1)
    {
        for (Int_t i_cent_pid = 1; i_cent_pid < vec_pid_cent_upper_dNdpt[i_particle].size(); i_cent_pid++)
        {
            if (sub_str == vec_pid_energy_dNdpt[i_particle][i_cent_pid])
            {
                ComboCentralityLabel_PID  = vec_pid_cent_lower_dNdpt[i_particle][i_cent_pid];
                ComboCentralityLabel_PID += "-";
                ComboCentralityLabel_PID += vec_pid_cent_upper_dNdpt[i_particle][i_cent_pid];
                ComboCentralityLabel_PID += " %";
                if ( i_number == 0 ) ComboCentrality_PID[i_type][i_particle]->AddEntry(ComboCentralityLabel_PID, i_cent_pid-1);
                vec_ComboCentralityLabel_PID.push_back(ComboCentralityLabel_PID.Copy());
                if (i_number>0)
                {
                    if (vec_ComboCentralityLabel_PID[i_number -1] != vec_ComboCentralityLabel_PID[i_number]) ComboCentrality_PID[i_type][i_particle]->AddEntry(ComboCentralityLabel_PID, i_cent_pid-1);
                }
                ComboCentralityLabel_PID.Clear();
                i_number ++;
            }
        }

    }
    return 1;
}
//______________________________________________________________________________
//______________________________________________________________________________
Int_t TBlastWaveGUI::DoSetSingleParticle(Int_t i_particle, Int_t i_type)
{
    // Called when energy and centrality is selected for a single particle.

    Int_t i_entry_energy = ComboEnergy_PID[i_type][i_particle]->GetSelected();
    if (i_entry_energy == -1) return 0;
    Int_t i_entry_centrality = ComboCentrality_PID[i_type][i_particle]->GetSelected();
    if (i_entry_centrality == -1) return 0;
    if (i_type == 0)
    {
        for ( Int_t i_entry = 0; i_entry < vec_tgae_id_v2.size(); i_entry++)
        {
            TString tgae_id_v2 = vec_tgae_id_v2[i_entry];
            Ssiz_t found = tgae_id_v2.First("_");
            tgae_id_v2.Replace(0, found+1, "");
            found = tgae_id_v2.First("_");
            tgae_id_v2.Replace(0, found+1, "");
            found = tgae_id_v2.First("_");
            TSubString sub_str =  tgae_id_v2(0,found);
            TString particle = sub_str;
            // cout <<"particle "  << particle <<label_full_pid_spectra[i_particle]  << endl;
            if ( particle  == label_full_pid_spectra[i_particle]) vec_tgae_id_v2.erase(vec_tgae_id_v2.begin()+i_entry);
        }
    }
    if ( i_type ==1)
    {
        for ( Int_t i_entry = 0; i_entry < vec_tgae_id_dNdpt.size(); i_entry++)
        {
            TString tgae_id_dNdpt = vec_tgae_id_dNdpt[i_entry];
            Ssiz_t found = tgae_id_dNdpt.First("_");
            tgae_id_dNdpt.Replace(0, found+1, "");
            found = tgae_id_dNdpt.First("_");
            tgae_id_dNdpt.Replace(0, found+1, "");
            found = tgae_id_dNdpt.First("_");
            TSubString sub_str =  tgae_id_dNdpt(0,found);
            TString particle = sub_str;
            // cout <<"particle "  << particle <<label_full_pid_spectra[i_particle]  << endl;
            if ( particle  == label_full_pid_spectra[i_particle]) vec_tgae_id_dNdpt.erase(vec_tgae_id_dNdpt.begin()+i_entry);
        }
    }


    ULong_t green;
    gClient->GetColorByName("green", green);

    if (i_type ==0)
    {
        fCheckBox_pid[i_particle]     ->ChangeBackground(green);
        fCheckBox_pid[i_particle]     ->SetState(kButtonDown);
        fCheckBox_pid_plot[i_particle]->ChangeBackground(green);
        fCheckBox_pid_plot[i_particle]->SetState(kButtonDown);
    }

    if (i_type ==1)
    {
        fCheckBox_pid_fit_dNdpt[i_particle] ->ChangeBackground(green);
        fCheckBox_pid_fit_dNdpt[i_particle] ->SetState(kButtonDown);
        fCheckBox_pid_plot_dNdpt[i_particle]->ChangeBackground(green);
        fCheckBox_pid_plot_dNdpt[i_particle]->SetState(kButtonDown);
    }

    TGTextLBEntry *cent_text_entry = (TGTextLBEntry *) ComboCentrality_PID[i_type][i_particle]->GetSelectedEntry();
    const char *selected_centrality = cent_text_entry->GetTitle();

    TString select_centrality = selected_centrality ;
    Ssiz_t found = select_centrality.First("-");
    TSubString sub_str = select_centrality(0,found);
    TString cent_lower = sub_str;

    select_centrality.Replace(0, found+1, "");
    found  = select_centrality.First(" ");
    sub_str = select_centrality(0,found);
    TString cent_upper = sub_str;

    TGTextLBEntry *energy_text_entry = (TGTextLBEntry *)ComboEnergy_PID[i_type][i_particle]->GetSelectedEntry();
    const char *selected_Energy = energy_text_entry->GetTitle();
    TString select_Energy = selected_Energy ;
    found = select_Energy.First(" ");
    TSubString sub_str_Energy = select_Energy(0,found);

    TString tgae_id;


    if (i_type ==0)
    {
        if (fCheckBox_pid_plot[i_particle]->IsDown() )   // pid plot v2
        {
            tgae_id = "v2_PID_";
            tgae_id += label_full_pid_spectra[i_particle].Data();
            tgae_id += "_E_";
            tgae_id += sub_str_Energy;
            tgae_id += "_C_";
            tgae_id += cent_lower.Data();
            tgae_id += "_";
            tgae_id += cent_upper.Data();
            cout << tgae_id << endl;
            vec_tgae_id_v2.push_back(tgae_id.Copy());
            tgae_id.Clear();
        }
        if (fCheckBox_pid[i_particle]->IsDown() )        // pid fit v2
        {
            tgae_id = "v2_PID_";
            tgae_id += label_full_pid_spectra[i_particle].Data();
            tgae_id += "_E_";
            tgae_id += sub_str_Energy;
            tgae_id += "_C_";
            tgae_id += cent_lower.Data();
            tgae_id += "_";
            tgae_id += cent_upper.Data();
            cout << tgae_id << endl;
            vec_tgae_id_v2_fit.push_back(tgae_id.Copy());
            tgae_id.Clear();
        }
    }
    if (i_type ==1)
    {
        if (fCheckBox_pid_plot_dNdpt[i_particle]->IsDown())    // pid plot dNdpt
        {
            tgae_id = "dNdpt_PID_";
            tgae_id += label_full_pid_spectra[i_particle].Data();
            tgae_id += "_E_";
            tgae_id += sub_str_Energy;
            tgae_id += "_C_";
            tgae_id += cent_lower.Data();
            tgae_id += "_";
            tgae_id += cent_upper.Data();
            cout << tgae_id << endl;
            vec_tgae_id_dNdpt.push_back(tgae_id.Copy());
            tgae_id.Clear();
        }

        if ( fCheckBox_pid_fit_dNdpt[i_particle]->IsDown())       // pid fit dNdpt
        {
            tgae_id = "dNdpt_PID_";
            tgae_id += label_full_pid_spectra[i_particle].Data();
            tgae_id += "_E_";
            tgae_id += sub_str_Energy;
            tgae_id += "_C_";
            tgae_id += cent_lower.Data();
            tgae_id += "_";
            tgae_id += cent_upper.Data();
            cout << tgae_id << endl;
            vec_tgae_id_dNdpt_fit.push_back(tgae_id.Copy());
            tgae_id.Clear();
        }
    }
   return 1;
}
//______________________________________________________________________________


//______________________________________________________________________________
void TBlastWaveGUI::PlotData()
{
    // Handle Plot data button.
    cout << "PlotData()" << endl;
    if(!c_1X1_v2)
    {
        c_1X1_v2 = new TCanvas("c_1X1_v2","c_1X1_v2",50,100,1400,900);
        c_1X1_v2->SetTopMargin(0.5);
        c_1X1_v2->SetBottomMargin(2);
        c_1X1_v2->SetRightMargin(0.5);
        c_1X1_v2->SetLeftMargin(1);
        c_1X1_v2->SetLogy(0);
    }

    c_1X1_v2->cd()->SetTicks(1,1);
    c_1X1_v2->cd()->SetGrid(0,0);
    c_1X1_v2->cd()->SetFillColor(10);
    c_1X1_v2->cd()->SetRightMargin(1);
    c_1X1_v2->cd()->SetTopMargin(1);

    h_dummy_plot_data->SetTitle("");
    h_dummy_plot_data->SetStats(0);
    h_dummy_plot_data->GetXaxis()->SetTitleOffset();
    h_dummy_plot_data->GetYaxis()->SetTitleOffset();
    h_dummy_plot_data->GetXaxis()->SetLabelOffset();
    h_dummy_plot_data->GetYaxis()->SetLabelOffset();
    h_dummy_plot_data->GetXaxis()->SetLabelSize();
    h_dummy_plot_data->GetYaxis()->SetLabelSize();
    h_dummy_plot_data->GetXaxis()->SetTitleSize();
    h_dummy_plot_data->GetYaxis()->SetTitleSize();
    h_dummy_plot_data->GetXaxis()->SetNdivisions(505,'N');
    h_dummy_plot_data->GetYaxis()->SetNdivisions(505,'N');
    h_dummy_plot_data->SetLineColor(10);
    h_dummy_plot_data->GetXaxis()->SetRangeUser(-0.03,10);
    h_dummy_plot_data->GetYaxis()->SetRangeUser(-0.1999,0.38);
    h_dummy_plot_data->DrawCopy();

    PlotLine(0.0,20,0,0,1,1,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle
    plotTopLegend((char*)"p_{T} (GeV/c)",0.475,0.03,0.05,1,0.0,42,1);
    plotTopLegend((char*)"#it{v}_{2}",0.03,0.5,0.05,1,90,42,1);


    
    if(!c_1X1_dNdpt)
    {
        c_1X1_dNdpt = new TCanvas("c_1X1_dNdpt","c_1X1_dNdpt",50,100,1400,900);
        c_1X1_dNdpt->SetTopMargin(0.5);
        c_1X1_dNdpt->SetBottomMargin(2);
        c_1X1_dNdpt->SetRightMargin(0.5);
        c_1X1_dNdpt->SetLeftMargin(1);
        c_1X1_dNdpt->SetLogy(0);
    }

    c_1X1_dNdpt->cd()->SetTicks(1,1);
    c_1X1_dNdpt->cd()->SetGrid(0,0);
    c_1X1_dNdpt->cd()->SetFillColor(10);
    c_1X1_dNdpt->cd()->SetRightMargin(1);
    c_1X1_dNdpt->cd()->SetTopMargin(1);

    h_dummy_plot_data->SetTitle("");
    h_dummy_plot_data->SetStats(0);
    h_dummy_plot_data->GetXaxis()->SetTitleOffset();
    h_dummy_plot_data->GetYaxis()->SetTitleOffset();
    h_dummy_plot_data->GetXaxis()->SetLabelOffset();
    h_dummy_plot_data->GetYaxis()->SetLabelOffset();
    h_dummy_plot_data->GetXaxis()->SetLabelSize();
    h_dummy_plot_data->GetYaxis()->SetLabelSize();
    h_dummy_plot_data->GetXaxis()->SetTitleSize();
    h_dummy_plot_data->GetYaxis()->SetTitleSize();
    h_dummy_plot_data->GetXaxis()->SetNdivisions(505,'N');
    h_dummy_plot_data->GetYaxis()->SetNdivisions(505,'N');
    h_dummy_plot_data->SetLineColor(10);
    h_dummy_plot_data->GetXaxis()->SetRangeUser(-0.03,6.5);
    h_dummy_plot_data->GetYaxis()->SetRangeUser(-0.1999,7);
    h_dummy_plot_data->DrawCopy();

    PlotLine(0.0,20,0,0,1,1,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle
    plotTopLegend((char*)"p_{T} (GeV/c)",0.475,0.03,0.05,1,0.0,42,1);
    plotTopLegend((char*)"dN/dp_{T}",0.03,0.5,0.05,1,90,42,1);

    c_1X1_v2 ->GetCanvas() ->cd();
    if(legend_v2_plot_data) delete legend_v2_plot_data;
    legend_v2_plot_data = new TLegend(0.9,0.5,0.8,0.82); // x1,y1,x2,y2
    legend_v2_plot_data->SetBorderSize(0);
    legend_v2_plot_data->SetFillColor(0);
    legend_v2_plot_data->SetTextSize(0.03);

    for (Int_t i_tgae_name = 0; i_tgae_name < (Int_t) vec_tgae_name_full.size(); i_tgae_name++)
    {
        Int_t i_mass = vec_index_pid[i_tgae_name];
        for (Int_t i_tgae_id = 0; i_tgae_id < (Int_t) vec_tgae_id_v2.size(); i_tgae_id++)
        {
            if ( vec_tgae_id_v2[i_tgae_id] == vec_tgae_name_full[i_tgae_name] && vec_error_type[i_tgae_name] =="stat")
            {
                if (fCheckBox_pid_plot[i_mass]->IsDown())
                {
                    vec_tgae[i_tgae_name]->SetMarkerSize(0.75);
                    vec_tgae[i_tgae_name]->SetMarkerStyle(20);
                    vec_tgae[i_tgae_name]->SetMarkerColor(arr_color_mass[i_mass]);
                    vec_tgae[i_tgae_name]->Draw("same PZ");
                    legend_v2_plot_data->AddEntry(vec_tgae[i_tgae_name],label_full_pid_spectra[i_mass],"p");
                    legend_v2_plot_data->Draw();
                }
            }

        }
    }

    c_1X1_dNdpt ->GetCanvas() ->cd();
    if(legend_dNdpt_plot_data) delete legend_dNdpt_plot_data;
    legend_dNdpt_plot_data = new TLegend(0.9,0.4,0.8,0.82); // x1,y1,x2,y2
    legend_dNdpt_plot_data->SetBorderSize(0);
    legend_dNdpt_plot_data->SetFillColor(0);
    legend_dNdpt_plot_data->SetTextSize(0.03);

    for (Int_t i_tgae_name = 0; i_tgae_name < (Int_t) vec_tgae_name_full.size(); i_tgae_name++)
    {
        Int_t i_mass = vec_index_pid[i_tgae_name];
        for (Int_t i_tgae_id = 0; i_tgae_id < (Int_t) vec_tgae_id_dNdpt.size(); i_tgae_id++)
        {
            if ( vec_tgae_id_dNdpt[i_tgae_id] == vec_tgae_name_full[i_tgae_name]&& vec_error_type[i_tgae_name] =="stat")
            {
                if (!fCheckBox_pid_plot_dNdpt[i_mass]->IsDown()) continue;
                vec_tgae[i_tgae_name]->SetMarkerSize(0.75);
                vec_tgae[i_tgae_name]->SetMarkerStyle(20);
                vec_tgae[i_tgae_name]->SetMarkerColor(arr_color_mass[i_mass]);
                vec_tgae[i_tgae_name]->Draw("same PZ");
                legend_dNdpt_plot_data->AddEntry(vec_tgae[i_tgae_name],label_full_pid_spectra[i_mass],"p");
                legend_dNdpt_plot_data->Draw();
            }
        }
    }
}

//______________________________________________________________________________
void TBlastWaveGUI::DoText(const char * /*text*/)
{
    //cout << "DoText()" << endl;
    // Handle text entry widgets.
    TGTextEntry *te = (TGTextEntry *) gTQSender;
    Int_t id = te->WidgetId();
    switch (id) {
    case HId1:
        fHslider1->SetPosition(atof(vec_TextBuffer[0]->GetString()),
                               fHslider1->GetMaxPosition());
        break;
    case HId2:
        fHslider1->SetPointerPosition(atof(fTbh2->GetString()));
        break;
    case HId3:
        fHslider1->SetPosition(fHslider1->GetMinPosition(),
                               atof(vec_TextBuffer[0]->GetString()));
        break;
    default:
        break;
    }
}



//______________________________________________________________________________
void TBlastWaveGUI::DoSlider()
{
    cout << "DoSlider started" << endl;

    Double_t number_entry = arr_NEntry_limits[0][0]->GetNumberEntry()->GetNumber();
    //printf("number_entry: %4.1f \n",number_entry);

    Int_t i_Temp   = vec_slider[0]->GetPosition();
    Int_t i_rho_0  = vec_slider[1]->GetPosition();
    Int_t i_rho_a  = vec_slider[2]->GetPosition();
    Int_t i_R_x    = vec_slider[3]->GetPosition();
    Int_t i_fboost = vec_slider[4]->GetPosition();

    Double_t Temp_val   = Temp_loop_start  + i_Temp*Delta_Temp;
    Double_t rho_0_val  = rho_0_loop_start + i_rho_0*Delta_rho_0;
    Double_t rho_a_val  = rho_a_loop_start + i_rho_a*Delta_rho_a;
    Double_t R_x_val    = arr_R_x[i_R_x];
    Double_t fboost_val = arr_f_boost[i_fboost];

    Double_t arr_param_val[5] = {Temp_val,rho_0_val,rho_a_val,R_x_val,fboost_val};

    for(Int_t i_mass = 0; i_mass < N_masses; i_mass++)
    {
        if(!tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]) continue;
        tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineWidth(3);
        tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineColor(arr_color_line_mass[i_mass]);
        tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineStyle(1);

        Double_t min_val_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();
        Double_t max_val_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();
        tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->GetXaxis()->SetRangeUser(min_val_pT,max_val_pT);

        h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineWidth(3);
        h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineColor(arr_color_line_mass[i_mass]);
        h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineStyle(1);
    }

    // Handle slider widgets.
    for(Int_t i_param = 0; i_param < 5; i_param++)
    {
        char buf[32];
        //sprintf(buf, "%d", vec_slider[i_param]->GetPosition());
        sprintf(buf,"%2.4f",arr_param_val[i_param]);
        vec_TextBuffer[i_param]->Clear();
        vec_TextBuffer[i_param]->AddText(0, buf);
        vec_TextEntry[i_param]->SetCursorPosition(vec_TextEntry[i_param]->GetCursorPosition());
        vec_TextEntry[i_param]->Deselect();
        gClient->NeedRedraw(vec_TextEntry[i_param]);
    }


    //-------------------------------------
    fCanvas ->GetCanvas() ->cd();
    h_dummy ->Draw();

    printf("Draw v2 \n");

    if(TL_line_base) delete TL_line_base;
    TL_line_base = PlotLine(0.0,15.0,0.0,0.0,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)

    for(Int_t i_mass = 0; i_mass < N_masses; i_mass++)
    {
        //cout << "Test A, i_mass: " << i_mass << endl;
        if(fCheckBox_pid_plot[i_mass]->GetState() == kButtonDown)
        {
            //printf("Draw i_mass: %d \n", i_mass);
            tgae_v2_vs_pT_mesons_data[i_mass] ->SetMarkerColor(arr_color_mass[i_mass]);
            tgae_v2_vs_pT_mesons_data[i_mass] ->SetLineColor(kGray);
            tgae_v2_vs_pT_mesons_data[i_mass] ->SetMarkerSize(0.8);
            tgae_v2_vs_pT_mesons_data[i_mass] ->SetMarkerStyle(20);
            tgae_v2_vs_pT_mesons_data[i_mass] ->Draw("same P");

            if(fCheckBox_sel[0]->GetState() == kButtonDown)
            {
                if(!tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]) continue;
                tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineWidth(4);
                tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineColor(arr_color_mass[i_mass]);
                tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->Draw("same L hist");

            }
        }
    }
    printf("v2 plotted \n");

    if(fCheckBox_sel[1]->GetState() == kButtonDown && flag_minimization_ana)
    {
        for(int i_mass = 0; i_mass < N_masses; ++i_mass)
        {
            if(fCheckBox_pid_plot[i_mass]->GetState() == kButtonDown)
            {
                tg_v2_BW_ana_pid[i_mass] -> Draw("same L");
            }
        }
    }

    if(leg_v2_vs_pT_A) delete leg_v2_vs_pT_A;
    leg_v2_vs_pT_A = new TLegend(0.76,0.63,0.86,0.82); // x1,y1,x2,y2
    leg_v2_vs_pT_A->SetBorderSize(0);
    leg_v2_vs_pT_A->SetFillColor(0);
    leg_v2_vs_pT_A->SetTextSize(0.045);
    leg_v2_vs_pT_A->AddEntry((TGraphAsymmErrors*)tgae_v2_vs_pT_mesons_data[0]->Clone(),"#pi","p");
    leg_v2_vs_pT_A->AddEntry((TGraphAsymmErrors*)tgae_v2_vs_pT_mesons_data[1]->Clone(),"K_{s}^{0}","p");
    leg_v2_vs_pT_A->AddEntry((TGraphAsymmErrors*)tgae_v2_vs_pT_mesons_data[2]->Clone(),"p","p");
    leg_v2_vs_pT_A->Draw();

    if(leg_v2_vs_pT_B) delete leg_v2_vs_pT_B;
    leg_v2_vs_pT_B = new TLegend(0.87,0.68,0.97,0.81); // x1,y1,x2,y2
    leg_v2_vs_pT_B->SetBorderSize(0);
    leg_v2_vs_pT_B->SetFillColor(0);
    leg_v2_vs_pT_B->SetTextSize(0.045);
    leg_v2_vs_pT_B->AddEntry((TGraphAsymmErrors*)tgae_v2_vs_pT_mesons_data[6]->Clone(),"J/#psi","p");
    leg_v2_vs_pT_B->AddEntry((TGraphAsymmErrors*)tgae_v2_vs_pT_mesons_data[7]->Clone(),"#varUpsilon","p");
    leg_v2_vs_pT_B->Draw();

    if(TLatex_legend_system[0]) delete TLatex_legend_system[0];
    TLatex_legend_system[0] = plotTopLegend((char*)"30-40%",0.73,0.83,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    if(TLatex_legend_system[1]) delete TLatex_legend_system[1];
    TLatex_legend_system[1] = plotTopLegend((char*)"5-60%",0.85,0.83,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    if(TLatex_legend_system[2]) delete TLatex_legend_system[2];
    TLatex_legend_system[2] = plotTopLegend((char*)"|y|<0.5",0.73,0.89,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    if(TLatex_legend_system[3]) delete TLatex_legend_system[3];
    TLatex_legend_system[3] = plotTopLegend((char*)"2.5<y<4",0.85,0.89,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    HistName = "#chi_{MC}^{2}/ndf = ";
    sprintf(NoP,"%4.2f",chi2_final_min);
    HistName += NoP;
    TLatex_chi2_MC = plotTopLegend((char*)HistName.Data(),0.18,0.86,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    HistName = "#chi_{Ana}^{2}/ndf = ";
    sprintf(NoP,"%4.2f",chi2_final_min_ana);
    HistName += NoP;
    if(TLatex_chi2_ana) delete TLatex_chi2_ana;
    TLatex_chi2_ana = plotTopLegend((char*)HistName.Data(),0.18,0.80,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1


    fCanvas->GetCanvas()->Modified();
    fCanvas->GetCanvas()->Update();
    //-------------------------------------


    //-------------------------------------
    Double_t x_range_dNdpT[N_masses] = {3.5,4.1,4.5,5.0,5.0,5.0,8.5,15.0,10.0};
    Double_t y_range_dNdpT[N_masses] = {2.2,1.4,1.1,1.0,1.0,1.0,0.45,0.35,0.5};

    for(Int_t i_mass = 0; i_mass < N_masses; i_mass++)
    {
        fCanvasB ->GetCanvas()->cd(i_mass + 1);
        h_dummy_dNdpT ->GetXaxis()->SetRangeUser(0.0,x_range_dNdpT[i_mass]);
        h_dummy_dNdpT ->GetYaxis()->SetRangeUser(-0.2,y_range_dNdpT[i_mass]);
        h_dummy_dNdpT ->DrawCopy("");

        tgae_dN_dpT_mesons_data[i_mass] ->SetMarkerColor(kBlack);
        tgae_dN_dpT_mesons_data[i_mass] ->SetMarkerStyle(20);
        tgae_dN_dpT_mesons_data[i_mass] ->SetMarkerSize(1.0);
        tgae_dN_dpT_mesons_data[i_mass] ->Draw("same P");

        if(!h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]) continue;
        h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineColor(kRed); // blast wave MC
        h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineWidth(5); // blast wave MC
        h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->DrawCopy("same hist L"); // blast wave MC
  
        if(flag_minimization_ana)
        {
            if(!tg_dNdpT_BW_ana_pid[i_mass]) continue;
            tg_dNdpT_BW_ana_pid[i_mass] -> SetLineColor(kAzure-2);
            tg_dNdpT_BW_ana_pid[i_mass] -> SetLineWidth(3);
            tg_dNdpT_BW_ana_pid[i_mass] -> SetLineStyle(1);
            if(!(fCheckBox_pid[i_mass]->GetState() == kButtonDown)) tg_dNdpT_BW_ana_pid[i_mass] -> SetLineStyle(9);
            if(fCheckBox_sel[1]->GetState() == kButtonDown)
            {
                tg_dNdpT_BW_ana_pid[i_mass] -> Draw("same L");
            }
        }
    }

    //------------------------------------------------------


    for(Int_t iPad = 1; iPad <= 8; iPad++)
    {
        fCanvasB ->GetCanvas()->cd(iPad);
        if(TLatex_legend_dNdpT[iPad-1]) delete TLatex_legend_dNdpT[iPad-1];
        TLatex_legend_dNdpT[iPad-1] = plotTopLegend((char*)label_pid_spectra[iPad-1].Data(),0.75,0.83,0.06,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    }

    fCanvasB->GetCanvas()->Modified();
    fCanvasB->GetCanvas()->Update();
    //-------------------------------------

}

//______________________________________________________________________________
void TBlastWaveGUI::StopMinimize()
{
    printf("Minimization stopped \n");
    flag_stop_minimize = 1;
}


//______________________________________________________________________________
void TBlastWaveGUI::DrawEllipse()
{
    cout << "DrawEllipse started" << endl;

    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_draw_ellipse_ana->ChangeBackground(green);

#if 0
    //------------------------------------------------------------
    printf("Do contour \n");
    if(!c_correlate)
    {
        c_correlate = new TCanvas("c_correlate","c_correlate",100,200,1400,600);
        c_correlate->SetFillColor(10);
        c_correlate->SetTopMargin(0.05);
        c_correlate->SetBottomMargin(0.22);
        c_correlate->SetRightMargin(0.05);
        c_correlate->SetLeftMargin(0.22);
        c_correlate->SetLogy(0);
    }


    unsigned int npoints = 100;
    double pntsxA[100];
    double pntsyA[100];
    resultC.Contour(0,1,npoints,pntsxA,pntsyA,0.683); // unsigned int ipar, unsigned int jpar, unsigned int &npoints, double *pntsx, double *pntsy, double confLevel=0.683
    gr1 = new TGraph(100,pntsxA,pntsyA);

    for(Int_t i_point = 0; i_point < npoints; i_point++)
    {
        printf("i_point: %d, x/y: {%4.3f, %4.3f} \n",i_point,pntsxA[i_point],pntsyA[i_point]);
    }

    /*
    double pntsxB[100];
    double pntsyB[100];
    result.Contour(0,1,100,pntsxB,pntsyB,0.955); // unsigned int ipar, unsigned int jpar, unsigned int &npoints, double *pntsx, double *pntsy, double confLevel=0.683
    gr2 = new TGraph(100,pntsxB,pntsyB);


    gr2->SetFillColor(42);
    gr2->Draw("lf");
    */
    c_correlate ->cd();
    gr1->SetFillColor(38);
    gr1->Draw("lf");

    printf("Contour done \n");
    //------------------------------------------------------------
#endif
}


//______________________________________________________________________________
void TBlastWaveGUI::DoMinimize_ana()
{
    cout << "DoMinimize_ana started" << endl;
    N_calls_BW_ana = 0;

    Int_t id_bw_hypersurface = fCombo ->GetSelected();
    printf("id_bw_hypersurface: %d \n",id_bw_hypersurface);

    fHProg2 ->Reset();

    Pixel_t yellow;
    gClient->GetColorByName("yellow", yellow);
    Button_minimize_ana->ChangeBackground(yellow);
    Button_take_params_MC_to_ana->ChangeBackground(yellow);

    //Particle Id list
    const Int_t idPi = 211;
    const Int_t idKaon = 321;
    const Int_t idProton = 2212;
    const Int_t idPhi = 333;
    const Int_t idXi = 3312;
    const Int_t idOmega = 3334;
    const Int_t idLambda = 3122;
    const Int_t idK0S = 310;
    const Int_t idD0 = 421;
    const Int_t idJPsi = 443;
    const Int_t idUpsilon = 553;
    const Int_t id_d = 0;
    const Int_t idHe3 =0;
    const Int_t id_t = 0;
    Int_t Partid[22] = {idPi,-idPi,idKaon,-idKaon,idProton,-idProton,idPhi,idXi,-idXi,idOmega,-idOmega,idLambda,-idLambda,idK0S,idD0,idJPsi,idUpsilon,0,0,0,0,0};
    Double_t Tch = 0.155;
    const Int_t index_mass_par = 0;
    const Int_t index_spintype_par = 1;
    /*
    TPythia8 tp;
    feeddown fd;
    fd.read_decay_histograms("./decay_histograms/hadron_decays_scan_01.root");
    fd.set_dndpt_and_v2_model(blastwave_dndpt_and_v2);
    fd.set_maximum_mother_mass(1.3);
    //-------------------------------------------------
    */
    // create and fill vectors with data for the fit
    vector <TGraphAsymmErrors*> tgae_v2_vs_pT_data;    // v2
    vector <TGraphAsymmErrors*> tgae_dN_dpT_data;      // dNdpt
    vector <Int_t> index_pid_v2_vs_pT_data;            // index of particles v2 (0-21)
    vector <Int_t> index_pid_dN_dpT_data;              // index of particles dNdpt (0-21)
    vector <TH1D*> dndpt_feed_hist;                    // dNdpt feeddown histogram
    vector <TH1D*> v2_feed_hist;                       // v2 feeddown histogram
    vector <Int_t> part_id;
    vector <Int_t> part_id_dndpt;
    vector <Int_t> part_id_v2;
    part_id.clear();
    part_id_dndpt.clear();
    part_id_v2.clear();
    index_pid_v2_vs_pT_data.clear();
    index_pid_dN_dpT_data.clear();
    tgae_v2_vs_pT_data.clear();
    tgae_dN_dpT_data.clear();


    for (Int_t i_tgae_name = 0; i_tgae_name < (Int_t) vec_tgae_name_full.size(); i_tgae_name++)
    {
        Int_t index_pid = vec_index_pid[i_tgae_name];
        if (!fCheckBox_pid[index_pid]->IsDown()) continue;  // only the particles with clicked PID fit v2 check boxes are considered
        for (Int_t i_tgae_id = 0; i_tgae_id < (Int_t) vec_tgae_id_v2_fit.size(); i_tgae_id++)
        {
            if ( vec_tgae_id_v2_fit[i_tgae_id] == vec_tgae_name_full[i_tgae_name] && vec_error_type[i_tgae_name] =="stat")
            {
                tgae_v2_vs_pT_data.push_back((TGraphAsymmErrors*)vec_tgae[i_tgae_name]->Clone());
                index_pid_v2_vs_pT_data.push_back(vec_index_pid[i_tgae_name]);
                part_id_v2.push_back(Partid[index_pid]);
            }
        }
    }

    for (Int_t i_tgae_name = 0; i_tgae_name < (Int_t) vec_tgae_name_full.size(); i_tgae_name++)
    {
        Int_t index_pid = vec_index_pid[i_tgae_name];
        if (!fCheckBox_pid_fit_dNdpt[index_pid]->IsDown()) continue;  // only the particles with clicked PID fit dNdpt check boxes are considered
        for (Int_t i_tgae_id = 0; i_tgae_id < (Int_t) vec_tgae_id_dNdpt_fit.size(); i_tgae_id++)
        {
            if ( vec_tgae_id_dNdpt_fit[i_tgae_id] == vec_tgae_name_full[i_tgae_name] && vec_error_type[i_tgae_name] =="stat")
            {
                tgae_dN_dpT_data.push_back((TGraphAsymmErrors*)vec_tgae[i_tgae_name]->Clone());
                index_pid_dN_dpT_data.push_back(vec_index_pid[i_tgae_name]);
                part_id_dndpt.push_back(Partid[index_pid]);
            }
        }
    }
    //--------------------------------------------------


    auto chi2Function = [&](const Double_t *par)
    {
        //minimisation function computing the sum of squares of residuals
        N_calls_BW_ana++;

        double chi2 = 0;

        const double T = par[0];        // fit parameter: Temp in GeV
        const double rho0 = par[1];     // fit parameter: transverse rapidity
        const double rho2 = par[2];     // fit parameter: azimuthal modulation of transverse rapidity
        const double RxOverRy = par[3]; // fit parameter: ratio of the radii Rx and Ry of the freeze-out ellipse in the transverse plane

        /*
        //feeddown calculations
        if(!(fCheckBoxFeedDown->GetState() == kButtonUp)){
            v2_feed_hist.clear();
            dndpt_feed_hist.clear();
            for(Int_t i_tgae = 0; i_tgae < (Int_t) part_id_dndpt.size(); i_tgae++)
            {
                Double_t model_pars[7] = {0., 0., T, rho0, rho2, RxOverRy, Tch};
                fd.set_model_parameters(model_pars,7);
                fd.calc_feeddown_hist(part_id_dndpt[i_tgae]);
                dndpt_feed_hist.push_back((TH1D*) fd.get_dndpt_total_hist()->Clone());
                v2_feed_hist.push_back((TH1D*) fd.get_v2_total_hist()->Clone());
            }
        }
        */
        // v2
        // loop over all graphs
        for(Int_t i_tgae = 0; i_tgae < (Int_t) tgae_v2_vs_pT_data.size(); i_tgae++)
        {
            Int_t i_mass_v2 = index_pid_v2_vs_pT_data[i_tgae];   // index to get the pT range (min_val_pT_v2, max_val_pT_v2)
            Double_t min_val_pT_v2;
            Double_t max_val_pT_v2;
            if(!(fCheckBox_pid[i_mass_v2]->GetState() == kButtonDown)) continue;
            if ( i_mass_v2 < 11)
            {
                min_val_pT_v2 = arr_NEntry_limits[0][i_mass_v2]->GetNumberEntry()->GetNumber();   // number entries for particles 0-10
                max_val_pT_v2 = arr_NEntry_limits[1][i_mass_v2]->GetNumberEntry()->GetNumber();
            }
            if ( i_mass_v2 >= 11)
            {
                min_val_pT_v2 = arr_NEntry_limits_A[0][i_mass_v2]->GetNumberEntry()->GetNumber();  // number entries for particle 11-21
                max_val_pT_v2 = arr_NEntry_limits_A[1][i_mass_v2]->GetNumberEntry()->GetNumber();
            }

            cout << "min_val_pT_v2: " << min_val_pT_v2 << "max_val_pT_v2: " << max_val_pT_v2  <<endl;
            individual_chi2_ana[i_mass_v2][0]   = 0.0; // v2
            N_individual_chi2_ana[i_mass_v2][0] = 0.0; // v2

            const double m = arr_quark_mass_meson[i_mass_v2];       // in GeV
            //cout << "mass: " << m <<endl;
            if(tg_v2_BW_ana_pid_min[i_mass_v2])    delete tg_v2_BW_ana_pid_min[i_mass_v2];

            tg_v2_BW_ana_pid_min[i_mass_v2]    = new TGraph();    // fit

            //--------------------------------------------------
            // v2 chi2
            Int_t i_point_ana = 0;
            // loop over all data points
            for(int i_point = 0; i_point < tgae_v2_vs_pT_data[i_tgae]->GetN(); ++i_point)
            {
                double v2_data = 0.0;
                double pt_data = 0.0;

                tgae_v2_vs_pT_data[i_tgae]->GetPoint(i_point,pt_data,v2_data);   // data
                double v2_err = tgae_v2_vs_pT_data[i_tgae]->GetErrorYhigh(i_point);
                cout << "i_point = " << i_point << ", pt_data = " << pt_data << ", v2_data = " << v2_data << " +/- " << v2_err << endl;

                if(pt_data > max_val_pT_v2) break;

                if(pt_data > min_val_pT_v2)
                {
                    double v2_BW = 0;
                    double inv_yield_BW = 0;

                    // blast wave parameters
                    const double pt_BW = pt_data;         // in GeV

                    /*
                    if(!(fCheckBoxFeedDown->GetState() == kButtonUp)){
                        TH1D *h_v2_fit = (TH1D*)v2_feed_hist[i_tgae];
                        Double_t v2_total = h_v2_fit->Interpolate(pt_BW);
                        v2_BW = v2_total;
                    }
                    */
                    if(!(fCheckBoxFeedDown->GetState() == kButtonDown) && id_bw_hypersurface == 1) bw_ana.calc_blastwave_yield_and_v2_fos1(pt_BW, m, T, rho0, rho2, RxOverRy, inv_yield_BW, v2_BW);
                    if(!(fCheckBoxFeedDown->GetState() == kButtonDown) && id_bw_hypersurface == 2) bw_ana.calc_blastwave_yield_and_v2_fos2(pt_BW, m, T, rho0, rho2, RxOverRy, inv_yield_BW, v2_BW);
                    if(!(fCheckBoxFeedDown->GetState() == kButtonDown) && id_bw_hypersurface == 3) bw_ana.calc_blastwave_yield_and_v2_fos3(pt_BW, m, T, rho0, rho2, RxOverRy, inv_yield_BW, v2_BW);
                    //blastwave_yield_and_v2(pt_BW, m, T, rho0, rho2, RxOverRy, inv_yield_BW, v2_BW);

                    tg_v2_BW_ana_pid_min[i_mass_v2] ->SetPoint(i_point_ana,pt_BW,v2_BW);

                    cout << "i_point = " << i_point <<  ", pt_BW = " << pt_BW << ", v2_BW = " << v2_BW<< ", chi2 = " << chi2 << ", T = " << T << ", rho0 = " << rho0 <<", rho2 = " << rho2 << ", RxOverRy = " << RxOverRy << endl;
                    double diff   = (v2_data - v2_BW)/v2_err;
                    //double diff_yield = (inv_yield_BW - inv_yield_data)/yield_err;
                    chi2 += diff*diff;
                    individual_chi2_ana[i_mass_v2][0]   += diff*diff;
                    N_individual_chi2_ana[i_mass_v2][0] += 1.0;
                    i_point_ana++;

                }
            }
        }

        //--------------------------------------------------


        //--------------------------------------------------
        // dNdpt
        //loop over all graphs in tgae_dN_dpT_data
        for(Int_t i_tgae = 0; i_tgae < (Int_t) tgae_dN_dpT_data.size(); i_tgae++)
        {
            Int_t i_mass = index_pid_dN_dpT_data[i_tgae];
            if(!(fCheckBox_pid_fit_dNdpt[i_mass]->GetState() == kButtonDown)) continue;

            // pT range
            Double_t min_val_pT;
            Double_t max_val_pT;
            if ( i_mass < 11)
            {
                min_val_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();   // number entries for particles 0-10
                max_val_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();
            }
            if ( i_mass >= 11)
            {
                min_val_pT = arr_NEntry_limits_A[0][i_mass]->GetNumberEntry()->GetNumber(); // number entries for particles 11-21
                max_val_pT = arr_NEntry_limits_A[1][i_mass]->GetNumberEntry()->GetNumber();
            }
            cout << "min_val_pT: " << min_val_pT<< "max_val_pT: " << max_val_pT  <<endl;


            individual_chi2_ana[i_mass][1]   = 0.0; // dN/dpT
            N_individual_chi2_ana[i_mass][1] = 0.0; // dN/dpT

            const double m = arr_quark_mass_meson[i_mass];       // in GeV

            if(tg_dNdpT_BW_ana_pid_min[i_mass]) delete tg_dNdpT_BW_ana_pid_min[i_mass];

            tg_dNdpT_BW_ana_pid_min[i_mass] = new TGraph();

            //--------------------------------------------------
            // dNdpT chi2

            Double_t integral_ana = 0.0;
            vector< vector<Double_t> > vec_data_BW;
            vec_data_BW.resize(4); // data, BW, err, pT

            // First the integral of BW needs to be calculated to do a shape comparison to the data
            // loop over all data points
            for(int i_point = 0; i_point < tgae_dN_dpT_data[i_tgae]->GetN(); ++i_point) 
            {
                double dNdpT_data = 0.0;
                double pt_data    = 0.0;

                tgae_dN_dpT_data[i_tgae]->GetPoint(i_point,pt_data,dNdpT_data);      // data dNdpt
                double dNdpT_err = tgae_dN_dpT_data[i_tgae]->GetErrorYhigh(i_point);

                Double_t x_err_low_data  = fabs(tgae_dN_dpT_data[i_tgae] ->GetErrorXlow(i_point));
                Double_t x_err_high_data = fabs(tgae_dN_dpT_data[i_tgae] ->GetErrorXhigh(i_point));
                Double_t bin_width = x_err_low_data + x_err_high_data;

                cout << "i_point = " << i_point << ", pt_data = " << pt_data << ", dNdpt_data = " << dNdpT_data << " +/- " << dNdpT_err << endl;

                double v2_BW = 0;
                double inv_yield_BW = 0;

                // blast wave parameters
                const double pt_BW = pt_data;         // in GeV

                /*
                if(!(fCheckBoxFeedDown->GetState() == kButtonUp)){
                    TH1D *h_dndpt_fit = (TH1D*)dndpt_feed_hist[i_tgae];
                    Double_t dndpt_total = h_dndpt_fit->Interpolate(pt_BW);;
                    inv_yield_BW = dndpt_total;
                }
                */
                if(!(fCheckBoxFeedDown->GetState() == kButtonDown) && id_bw_hypersurface == 1) bw_ana.calc_blastwave_yield_and_v2_fos1(pt_BW, m, T, rho0, rho2, RxOverRy, inv_yield_BW, v2_BW); // invariant yield: (1/pT) (dN/dpT)
                if(!(fCheckBoxFeedDown->GetState() == kButtonDown) && id_bw_hypersurface == 2) bw_ana.calc_blastwave_yield_and_v2_fos2(pt_BW, m, T, rho0, rho2, RxOverRy, inv_yield_BW, v2_BW); // invariant yield: (1/pT) (dN/dpT)
                if(!(fCheckBoxFeedDown->GetState() == kButtonDown) && id_bw_hypersurface == 3) bw_ana.calc_blastwave_yield_and_v2_fos3(pt_BW, m, T, rho0, rho2, RxOverRy, inv_yield_BW, v2_BW); // invariant yield: (1/pT) (dN/dpT)
                vec_data_BW[0].push_back(dNdpT_data);
                vec_data_BW[1].push_back(inv_yield_BW);
                vec_data_BW[2].push_back(dNdpT_err);
                vec_data_BW[3].push_back(pt_BW);

                integral_ana += inv_yield_BW*bin_width;
                cout << "i_point = " << i_point << ", pt_BW = " << pt_BW << ", v2_BW = " << v2_BW << ", T = " << T <<", rho0 = " << rho0 <<", rho2 = " << rho2 <<", RxOverRy = " << RxOverRy <<", inv_yield_BW = " << inv_yield_BW << endl;
                //double diff   = (dNdpT_data - inv_yield_BW)/dNdpT_err;
                //double diff_yield = (inv_yield_BW - inv_yield_data)/yield_err;
                //chi2 += diff*diff;
            }

            Double_t scale_factor_ana = get_norm_scaling_factor_calc(vec_data_BW,min_val_pT,max_val_pT);

            // Calculate chi2 for dNdpT
            if(scale_factor_ana > 0.0)
            {
                Int_t i_point_ana = 0;
                for(Int_t i_point = 0; i_point < (Int_t)vec_data_BW[0].size(); i_point++)
                {
                    Double_t pt_data = vec_data_BW[3][i_point];
                    if(pt_data > max_val_pT) break;

                    if(pt_data > min_val_pT)
                    {
                        // Normalize BW to integral within range of data
                        double diff   = (vec_data_BW[0][i_point] - (vec_data_BW[1][i_point]*scale_factor_ana))/vec_data_BW[2][i_point];
                        chi2 += diff*diff;
                        individual_chi2_ana[i_mass][1]   += diff*diff;
                        N_individual_chi2_ana[i_mass][1] += 1.0;
                        tg_dNdpT_BW_ana_pid_min[i_mass] ->SetPoint(i_point_ana,pt_data,vec_data_BW[1][i_point]*scale_factor_ana);
                        i_point_ana++;
                        cout << "i_point = " << i_point << ", chi2= " << chi2 << endl;
                    }
                }
            }

            //integration_range_pid[i_mass][0] = x_val_data_first - x_err_low_data;
            //integration_range_pid[i_mass][1] = x_val_data_last  + x_err_high_data;
            //--------------------------------------------------
        }

        if(N_calls_BW_ana % 5 == 0) fHProg2->Increment(5);
        printf("N_calls_BW_ana: %d \n",N_calls_BW_ana);
        //printf("fraction total: %4.2f%%, fraction added: %4.2f%% \n",fraction_total*100.0,fraction_use*100.0);

        // v2 fit plot
        c_1X1_v2->GetCanvas()->cd();
        HistName = "#chi_{Ana}^{2}/ndf = ";
        sprintf(NoP,"%4.2f",chi2);
        HistName += NoP;
        //printf("chi2 (ana): %4.3f \n",chi2);
        if(TLatex_chi2_ana) delete TLatex_chi2_ana;
        TLatex_chi2_ana = plotTopLegend((char*)HistName.Data(),0.18,0.80,0.03,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        for(int i_mass = 0; i_mass < N_masses_all; ++i_mass)
        {
            if(!tg_v2_BW_ana_pid_min[i_mass]) continue;
            tg_v2_BW_ana_pid_min[i_mass] ->SetLineColor(arr_color_mass[i_mass]);
            tg_v2_BW_ana_pid_min[i_mass] ->SetLineWidth(3);
            tg_v2_BW_ana_pid_min[i_mass] ->SetLineStyle(9);
            tg_v2_BW_ana_pid_min[i_mass] ->Draw("same L");
        }

        c_1X1_v2->GetCanvas()->Modified();
        c_1X1_v2->GetCanvas()->Update();

        // dNdpt fit plot
        c_1X1_dNdpt ->GetCanvas()->cd();
        for(int i_mass = 0; i_mass < N_masses_all; ++i_mass)
        {
            if(!tg_dNdpT_BW_ana_pid_min[i_mass]) continue;
            //fCanvasB ->GetCanvas()->cd(i_mass+1);
            tg_dNdpT_BW_ana_pid_min[i_mass] ->SetLineColor(arr_color_mass[i_mass]);
            tg_dNdpT_BW_ana_pid_min[i_mass] ->SetLineWidth(3);
            tg_dNdpT_BW_ana_pid_min[i_mass] ->SetLineStyle(1);
            tg_dNdpT_BW_ana_pid_min[i_mass] ->Draw("same L");
        }

        c_1X1_dNdpt ->GetCanvas()->Modified();
        c_1X1_dNdpt ->GetCanvas()->Update();

        gSystem->Sleep(100);
        gSystem->ProcessEvents();

        if(chi2 < chi2_final_min_ana)
        {
            chi2_final_min_ana = chi2;
            for(int i_mass = 0; i_mass < N_masses_all; ++i_mass)
            {
                if(N_individual_chi2_ana[i_mass][0] > 0.0)
                {
                    best_individual_chi2_per_point_ana[i_mass][0] = individual_chi2_ana[i_mass][0]/N_individual_chi2_ana[i_mass][0];
                }
                else best_individual_chi2_per_point_ana[i_mass][0] = -1.0;
                if(N_individual_chi2_ana[i_mass][1] > 0.0)
                {
                    best_individual_chi2_per_point_ana[i_mass][1] = individual_chi2_ana[i_mass][1]/N_individual_chi2_ana[i_mass][1];
                }
                else best_individual_chi2_per_point_ana[i_mass][1] = -1.0;
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
    //fitter.Config().ParSettings(0).Fix(); // T
    fitter.Config().ParSettings(0).SetLimits(0.05,2.5); // T
    fitter.Config().ParSettings(1).SetLimits(0.01,2.0); // rho0
    fitter.Config().ParSettings(2).SetLimits(0.01,1.0); // rho2
    fitter.Config().ParSettings(3).SetLimits(0.01,2.0); // RxOverRy

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

    const ROOT::Fit::FitResult &resultB = fitter.Result();
    ROOT::Fit::FitResult result(resultB);
    resultC = result;
    result.Print(std::cout);
    const double *fitpar = result.GetParams();
    double T_BW        = fitpar[0];
    double Rho0_BW     = fitpar[1];
    double Rho2_BW     = fitpar[2];
    double RxOverRy_BW = fitpar[3];

    const double *fitpar_err = result.GetErrors();
    double T_BW_err        = fitpar_err[0];
    double Rho0_BW_err     = fitpar_err[1];
    double Rho2_BW_err     = fitpar_err[2];
    double RxOverRy_BW_err = fitpar_err[3];

    T_BW_fit_ana        = T_BW;
    Rho0_BW_fit_ana     = Rho0_BW;
    Rho2_BW_fit_ana     = Rho2_BW;
    RxOverRy_BW_fit_ana = RxOverRy_BW;
    Double_t mean_beta = TMath::TanH(Rho0_BW_fit_ana);


    double Chi2_ana = result.Chi2();
    double NDF_ana  = result.Ndf();
    cout << "T_BW = " << T_BW << ", Rho0_BW = " << Rho0_BW << ", Rho2_BW = " << Rho2_BW << ", RxOverRy_BW = " << RxOverRy_BW << ", mean_beta: " << mean_beta << endl;

    TString par_names[4] = {"T","rho0","rho2","Rx"};

    printf("%10s %10s %10s %10s %10s \n"," ","T","rho0","rho2","Rx/Ry");

    for(Int_t i_parA = 0; i_parA < 4; i_parA++)
    {
        printf("%10s",par_names[i_parA].Data());
        for(Int_t i_parB = 0; i_parB < 4; i_parB++)
        {
            Double_t cov_value =  result.CovMatrix(i_parA,i_parB);
            printf("%10.5f",cov_value/(fitpar_err[i_parA]*fitpar_err[i_parB]));
            if(i_parB == 3) printf(" \n");
            //printf("%s-%s: %4.7f \n",par_names[i_parA].Data(),par_names[i_parB].Data(),cov_value);
        }
    }

    //Plot_curves_ana(T_BW,Rho0_BW,Rho2_BW,RxOverRy_BW);

    //-------------------------------------------------------------
    // save fit parameters and pT ranges (see function WriteParams)
    fit_params[0] = T_BW;
    fit_params[1] = Rho0_BW;
    fit_params[2] = Rho2_BW;
    fit_params[3] = RxOverRy_BW;

    Double_t min_val_pT, max_val_pT;
    for ( Int_t i_mass =0; i_mass < N_masses_all; i_mass++)
    {
        if ( i_mass < 11)
        {
            min_val_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();
            max_val_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();
        }
        if ( i_mass >= 11)
        {
            min_val_pT = arr_NEntry_limits_A[0][i_mass]->GetNumberEntry()->GetNumber();
            max_val_pT = arr_NEntry_limits_A[1][i_mass]->GetNumberEntry()->GetNumber();
        }
        vec_min_val_pT.push_back(min_val_pT);
        vec_max_val_pT.push_back(max_val_pT);
    }


    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_minimize_ana->ChangeBackground(green);


#if 0
    //------------------------------------------------------------
    printf("Do contour \n");
    if(!c_correlate)
    {
        c_correlate = new TCanvas("c_correlate","c_correlate",100,200,1400,600);
        c_correlate->SetFillColor(10);
        c_correlate->SetTopMargin(0.05);
        c_correlate->SetBottomMargin(0.22);
        c_correlate->SetRightMargin(0.05);
        c_correlate->SetLeftMargin(0.22);
        c_correlate->SetLogy(0);
    }


    unsigned int npoints = 50;
    double pntsxA[50];
    double pntsyA[50];
    result.Contour(2,3,npoints,pntsxA,pntsyA,0.683); // unsigned int ipar, unsigned int jpar, unsigned int &npoints, double *pntsx, double *pntsy, double confLevel=0.683
    gr1 = new TGraph(npoints,pntsxA,pntsyA);

    for(Int_t i_point = 0; i_point < npoints; i_point++)
    {
        printf("i_point: %d, x/y: {%4.3f, %4.3f} \n",i_point,pntsxA[i_point],pntsyA[i_point]);
    }


    double pntsxB[50];
    double pntsyB[50];
    result.Contour(2,3,npoints,pntsxB,pntsyB,0.8); // unsigned int ipar, unsigned int jpar, unsigned int &npoints, double *pntsx, double *pntsy, double confLevel=0.683
    gr2 = new TGraph(npoints,pntsxB,pntsyB);


    gr2->SetFillColor(42);
    gr2->Draw("lf");

    c_correlate ->cd();
    gr1->SetFillColor(38);
    gr1->Draw("lf");

    c_correlate->GetCanvas()->Modified();
    c_correlate->GetCanvas()->Update();

    printf("Contour done \n");
    //------------------------------------------------------------
#endif


}



//______________________________________________________________________________
void TBlastWaveGUI::DoMinimize()
{
    cout << "DoMinimize started" << endl;

    Pixel_t yellow;
    gClient->GetColorByName("yellow", yellow);
    Button_minimize->ChangeBackground(yellow);

    Button_take_params_MC_to_ana->ChangeBackground(yellow);

    flag_stop_minimize = 0;


    fHProg1 ->Reset();

    // new
    // Set pT min max values to those entered
    for(Int_t i_min_max = 0; i_min_max < 2; i_min_max++)
    {
        for(Int_t i_pid = 0; i_pid < N_masses; i_pid++)
        {
            Double_t number_entry = arr_NEntry_limits[i_min_max][i_pid]->GetNumberEntry()->GetNumber();
            min_max_pT_range_pid[i_min_max][i_pid] = number_entry;
        }
    }

    Int_t i_Temp_min   = 0;
    Int_t i_rho_0_min  = 0;
    Int_t i_rho_a_min  = 0;
    Int_t i_R_x_min    = 0;
    Int_t i_fboost_min = 0;

    Int_t i_Temp_dNdpT_min   = 0;
    Int_t i_rho_0_dNdpT_min  = 0;
    Int_t i_rho_a_dNdpT_min  = 0;
    Int_t i_R_x_dNdpT_min    = 0;
    Int_t i_fboost_dNdpT_min = 0;

    chi2_min                       = 10000000.0;
    Double_t chi2_tot              = 0;
    Double_t chi2_tot_piKp         = 0;
    chi2_ndf_dNdpT_min             = 10000000.0;
    Double_t chi2_ndf_dNdpT_min_v2 = 0;
    Double_t chi2_final            = 0;
    Double_t chi2_final_piKp       = 0;
    chi2_final_min                 = 10000000.0;
    chi2_final_min_piKp            = 10000000.0;
    Double_t chi2_ndf_dNdpT        = 0;
    Double_t chi2_ndf_dNdpT_piKp   = 0;

    Double_t N_total_params = TMath::Power(9,5);

    Double_t N_params_use = 0;
    Int_t    N_params_total_use = 0;
    Double_t fraction_progress_bar_update = 0.005;

    Int_t n_arr;
    Int_t n_arr_dNdpT;

    for (Int_t i_R_x = 0; i_R_x < 9; i_R_x++)
    {
        if(flag_stop_minimize) break;
        for (Int_t i_fboost = 0; i_fboost < 9; i_fboost++)
        {
            for (Int_t i_Temp = 0; i_Temp < 9; i_Temp++)
            {
                if(flag_stop_minimize) break;
                //printf("i_Temp: %d \n",i_Temp);
                for (Int_t i_rho_0 = 0; i_rho_0 < 9; i_rho_0++)
                {
                    if(flag_stop_minimize) break;
                    //printf("i_rho_0: %d \n",i_rho_0);
                    for (Int_t i_rho_a = 0; i_rho_a < 9; i_rho_a++)
                    {
                        if(flag_stop_minimize) break;
                        //printf("i_rho_a: %d \n",i_rho_a);

                        if(flag_stop_minimize) break;

                        Int_t    nop_tot        = 0; // number of points - sum
                        Int_t    nop_tot_piKp   = 0;
                        Int_t    nop_dNdpT      = 0;
                        Int_t    nop_dNdpT_piKp = 0;

                        chi2_final     = 0;
                        chi2_tot       = 0;
                        chi2_ndf_dNdpT = 0;


                        for(Int_t i_mass = 0; i_mass < N_masses; i_mass++) //

                            // dNdpT chi2 separately: pi good, K good, p very good,
                            //phi ~159, Omega very bad ~100000, D0 ~2000, J/Psi ~159, Upsilon (in fact in function_BW it's D again) ~697

                        {
                            if(!(fCheckBox_pid[i_mass]->GetState() == kButtonDown)) continue;

                            n_arr       = 0;
                            n_arr_dNdpT = 0;
                            n_arr       = tgae_v2_vs_pT_mesons_data[i_mass] ->GetN();
                            n_arr_dNdpT = tgae_dN_dpT_mesons_data[i_mass]   ->GetN();


                            Double_t x_pid;
                            Double_t v2_pid;
                            Double_t v2_err_pid;
                            Double_t v2_bw_pid;

                            if((fCheckBox_v2_dNdpT[0]->GetState() == kButtonDown))
                            {
                                for(Int_t i_pT = 0; i_pT < n_arr; i_pT++) // pT loop
                                {
                                    //get v2_pid, v2_err_pid for different particles

                                    tgae_v2_vs_pT_mesons_data[i_mass]              ->GetPoint(i_pT,x_pid,v2_pid);
                                    v2_err_pid = tgae_v2_vs_pT_mesons_data[i_mass] ->GetErrorY(i_pT);

                                    if(x_pid >= min_max_pT_range_pid[0][i_mass]
                                       && x_pid <= min_max_pT_range_pid[1][i_mass]
                                       && v2_err_pid != 0) // calculate only within cetrain pT range; one loop for everybody
                                    {
                                        v2_bw_pid           = tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->GetBinContent(tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]->FindBin(x_pid));
                                        chi2_tot            += ((v2_pid-v2_bw_pid)*(v2_pid-v2_bw_pid))/(v2_err_pid*v2_err_pid);
                                        nop_tot             += 1;

                                        // if (i_mass < 3) // pi K p only
                                        // {
                                        //     chi2_tot_piKp           += ((v2_pid-v2_bw_pid)*(v2_pid-v2_bw_pid))/(v2_err_pid*v2_err_pid);
                                        //     nop_tot_piKp            += 1;
                                        //}
                                    }

                                }
                            }

                            x_pid = 0;
                            Double_t dNdpT_pid;
                            Double_t dNdpT_err_pid;
                            Double_t dNdpT_bw_pid;

                            //if(i_mass != 7) // Upsilon spectra do not exist
                            {
                                if(!h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]) continue;
                                if((fCheckBox_v2_dNdpT[1]->GetState() == kButtonDown))
                                {
                                    for(Int_t i_pT = 0; i_pT < n_arr_dNdpT; i_pT++) // pT loop for dNdpT
                                    {
                                        //get dNdpT + errors

                                        tgae_dN_dpT_mesons_data[i_mass]                 ->GetPoint(i_pT,x_pid,dNdpT_pid);
                                        dNdpT_err_pid = tgae_dN_dpT_mesons_data[i_mass] ->GetErrorY(i_pT);

                                        if(x_pid >= min_max_pT_range_pid[0][i_mass]
                                           && x_pid <= min_max_pT_range_pid[1][i_mass]
                                           && dNdpT_err_pid != 0)  //&& (dNdpT_pid-dNdpT_bw_pid) <  100.0
                                        {
                                            dNdpT_bw_pid    = h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->GetBinContent(h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]->FindBin(x_pid));
                                            chi2_ndf_dNdpT += ((dNdpT_pid-dNdpT_bw_pid)*(dNdpT_pid-dNdpT_bw_pid))/(dNdpT_err_pid*dNdpT_err_pid);
                                            nop_dNdpT      += 1;

                                        }


                                    }
                                }
                            }


                        } // end of mass loop

                        chi2_final      = chi2_tot + chi2_ndf_dNdpT;
                        //chi2_final_piKp = chi2_tot_piKp + chi2_ndf_dNdpT_piKp;

                        if (nop_tot + nop_dNdpT > 5) //sum of chi2
                        {
                            chi2_final = chi2_final/(nop_tot + nop_dNdpT - 5);
                        }

                        // if (nop_tot_piKp + nop_dNdpT_piKp > 5) //sum of chi2
                        // {
                        //     chi2_final_piKp = chi2_final_piKp/(nop_tot_piKp + nop_dNdpT_piKp - 5);
                        // }


                        if(chi2_final < chi2_final_min)
                        {
                            chi2_final_min     = chi2_final;
                            i_Temp_min         = i_Temp;
                            i_rho_0_min        = i_rho_0;
                            i_rho_a_min        = i_rho_a;
                            i_R_x_min          = i_R_x;
                            i_fboost_min       = i_fboost;

                            printf("chi2_min: %4.3f \n",chi2_final_min);

                            vec_slider[0]->SetPosition(i_Temp_min);
                            vec_slider[1]->SetPosition(i_rho_0_min);
                            vec_slider[2]->SetPosition(i_rho_a_min);
                            vec_slider[3]->SetPosition(i_R_x_min);
                            vec_slider[4]->SetPosition(i_fboost_min);
                            gSystem->ProcessEvents();
                            DoSlider();
                            gSystem->Sleep(100);
                            gSystem->ProcessEvents();
                        }


                        N_params_use += 1.0;
                        N_params_total_use++;
                        Double_t fraction_use = (N_params_use/N_total_params);
                        Double_t fraction_total = (Double_t)N_params_total_use/N_total_params;
                        if(fraction_use >= fraction_progress_bar_update)
                        {
                            N_params_use = 0.0;
                            fHProg1->Increment(fraction_use*100.0);
                            printf("fraction total: %4.2f%%, fraction added: %4.2f%% \n",fraction_total*100.0,fraction_use*100.0);
                            gSystem->Sleep(100);
                            gSystem->ProcessEvents();
                            //cout << "End of sleep" << endl;
                        }
                    }
                }
            }
        }
    }

    if(flag_stop_minimize) cout << "Minimization was externally stopped" << endl;

    Double_t Temp_val_min   = Temp_loop_start  + i_Temp_min*Delta_Temp;
    Double_t rho_0_val_min  = rho_0_loop_start + i_rho_0_min*Delta_rho_0;
    Double_t rho_a_val_min  = rho_a_loop_start + i_rho_a_min*Delta_rho_a;
    Double_t R_x_val_min    = arr_R_x[i_R_x_min];
    Double_t fboost_val_min = arr_f_boost[i_fboost_min];

    Double_t arr_param_val_min[5] = {Temp_val_min,rho_0_val_min,rho_a_val_min,R_x_val_min,fboost_val_min};

    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_minimize->ChangeBackground(green);

    cout << "Hello Dave. Your pi K p minimum Chi2 = " << chi2_final_min << endl;
    cout << "Corresponding parameters are:" << endl;
    cout << "Temperature: " << arr_param_val_min[0] << endl;
    cout << "rho_0: " << arr_param_val_min[1] << endl;
    cout << "rho_a: " << arr_param_val_min[2] << endl;
    cout << "R_x: " << arr_param_val_min[3] << endl;
    cout << "fboost: " << arr_param_val_min[4] << endl;
}


//______________________________________________________________________________
void TBlastWaveGUI::TakeParamsFromSet()
{
    printf("TBlastWaveGUI::TakeParamsFromSet() \n");
    //flag_minimization_ana = 1;

    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_take_params_Set_to_ana->ChangeBackground(green);

    Pixel_t yellow;
    gClient->GetColorByName("yellow", yellow);
    Button_take_params_MC_to_ana->ChangeBackground(yellow);

    Double_t Temp = arr_NEntry_ana_params[0]->GetNumberEntry()->GetNumber();
    Double_t rho0 = arr_NEntry_ana_params[1]->GetNumberEntry()->GetNumber();
    Double_t rhoa = arr_NEntry_ana_params[2]->GetNumberEntry()->GetNumber();
    Double_t Rx   = arr_NEntry_ana_params[3]->GetNumberEntry()->GetNumber();

    T_BW_fit_ana        = Temp;
    Rho0_BW_fit_ana     = rho0;
    Rho2_BW_fit_ana     = rhoa;
    RxOverRy_BW_fit_ana = Rx;

    Plot_curves_ana(Temp,rho0,rhoa,Rx);
}


//______________________________________________________________________________
void TBlastWaveGUI::TakeParamsFromMC()
{
    printf("TBlastWaveGUI::TakeParamsFromMC() \n");
    flag_minimization_ana = 1;

    Int_t i_Temp   = vec_slider[0]->GetPosition();
    Int_t i_rho_0  = vec_slider[1]->GetPosition();
    Int_t i_rho_a  = vec_slider[2]->GetPosition();
    Int_t i_R_x    = vec_slider[3]->GetPosition();
    Int_t i_fboost = vec_slider[4]->GetPosition();

    Double_t Temp_val   = Temp_loop_start  + i_Temp*Delta_Temp;
    Double_t rho_0_val  = rho_0_loop_start + i_rho_0*Delta_rho_0;
    Double_t rho_a_val  = rho_a_loop_start + i_rho_a*Delta_rho_a;
    Double_t R_x_val    = arr_R_x[i_R_x];
    Double_t fboost_val = arr_f_boost[i_fboost];

    T_BW_fit_ana        = Temp_val;
    Rho0_BW_fit_ana     = rho_0_val;
    Rho2_BW_fit_ana     = rho_a_val;
    RxOverRy_BW_fit_ana = R_x_val;

    Plot_curves_ana(Temp_val,rho_0_val,rho_a_val,R_x_val);


    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_take_params_MC_to_ana->ChangeBackground(green);

    Pixel_t yellow;
    gClient->GetColorByName("yellow", yellow);
    Button_take_params_Set_to_ana->ChangeBackground(yellow);

    fCanvas->GetCanvas()->Modified();
    fCanvas->GetCanvas()->Update();
}
//______________________________________________________________________________



//______________________________________________________________________________
void TBlastWaveGUI::Plot_curves_ana(Double_t T_BW,Double_t  Rho0_BW,Double_t  Rho2_BW,Double_t  RxOverRy_BW)
{
    printf("TBlastWaveGUI::Plot_curves_ana() \n");
    flag_minimization_ana = 1;

    Int_t id_bw_hypersurface = fCombo ->GetSelected();
    printf("id_bw_hypersurface: %d \n",id_bw_hypersurface);

    //--------------------------------------------------------------------
    Double_t min_max_pT_range_pid_plot_ana[2][N_masses_all];
    min_max_pT_range_pid_plot_ana[0][0] = 0.1; // pi+
    min_max_pT_range_pid_plot_ana[1][0] = 3.5;
    min_max_pT_range_pid_plot_ana[0][1] = 0.1; // pi-
    min_max_pT_range_pid_plot_ana[1][1] = 3.5;
    min_max_pT_range_pid_plot_ana[0][2] = 0.1; // K+
    min_max_pT_range_pid_plot_ana[1][2] = 4.0;
    min_max_pT_range_pid_plot_ana[0][3] = 0.1; // K-
    min_max_pT_range_pid_plot_ana[1][3] = 4.0;
    min_max_pT_range_pid_plot_ana[0][4] = 0.1; // p
    min_max_pT_range_pid_plot_ana[1][4] = 3.0;
    min_max_pT_range_pid_plot_ana[0][5] = 0.1; // pbar
    min_max_pT_range_pid_plot_ana[1][5] = 3.0;
    min_max_pT_range_pid_plot_ana[0][6] = 0.1; // phi
    min_max_pT_range_pid_plot_ana[1][6] = 5.0;
    min_max_pT_range_pid_plot_ana[0][7] = 0.1; // xi-
    min_max_pT_range_pid_plot_ana[1][7] = 5.0;
    min_max_pT_range_pid_plot_ana[0][8] = 0.1; // xibar+
    min_max_pT_range_pid_plot_ana[1][8] = 5.0;
    min_max_pT_range_pid_plot_ana[0][9] = 0.1; // Omega-
    min_max_pT_range_pid_plot_ana[1][9] = 6.0;
    min_max_pT_range_pid_plot_ana[0][10] = 0.1; // Omegabar+
    min_max_pT_range_pid_plot_ana[1][10] = 6.0;
    min_max_pT_range_pid_plot_ana[0][11] = 0.1; // Lambda
    min_max_pT_range_pid_plot_ana[1][11] = 6.0;
    min_max_pT_range_pid_plot_ana[0][12] = 0.1; //Lambdabar
    min_max_pT_range_pid_plot_ana[1][12] = 6.0;
    min_max_pT_range_pid_plot_ana[0][13] = 0.1; // K0S
    min_max_pT_range_pid_plot_ana[1][13] = 4.0;
    min_max_pT_range_pid_plot_ana[0][14] = 0.1; // D0
    min_max_pT_range_pid_plot_ana[1][14] = 12.5;
    min_max_pT_range_pid_plot_ana[0][15] = 0.1; // J/Psi
    min_max_pT_range_pid_plot_ana[1][15] = 15.0;
    min_max_pT_range_pid_plot_ana[0][16] = 0.1; // Upsilon
    min_max_pT_range_pid_plot_ana[1][16] = 15.0;
    min_max_pT_range_pid_plot_ana[0][17] = 0.1; // d
    min_max_pT_range_pid_plot_ana[1][17] = 15.0;
    min_max_pT_range_pid_plot_ana[0][18] = 0.1; // dbar
    min_max_pT_range_pid_plot_ana[1][18] = 5.0;
    min_max_pT_range_pid_plot_ana[0][19] = 0.1; // He3
    min_max_pT_range_pid_plot_ana[1][19] = 5.0;
    min_max_pT_range_pid_plot_ana[0][20] = 0.1; // He3bar
    min_max_pT_range_pid_plot_ana[1][20] = 5.0;
    min_max_pT_range_pid_plot_ana[0][21] = 0.1; // t
    min_max_pT_range_pid_plot_ana[1][21] = 5.0;
    //--------------------------------------------------------------------


    //------------------------------------------------------
    // v2 plots
    c_1X1_v2 ->GetCanvas() ->cd();

    const Int_t N_points_BW_ana = 35;
    for(int i_mass = 0; i_mass < N_masses_all; ++i_mass)
    {
        if(!(fCheckBox_pid[i_mass]->GetState() == kButtonDown)) continue;
        Double_t delta_pT = (Double_t)((min_max_pT_range_pid_plot_ana[1][i_mass] - min_max_pT_range_pid_plot_ana[0][i_mass])/(Double_t)(N_points_BW_ana-1));
        //Double_t min_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();
        //Double_t max_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();


        if(tg_v2_BW_ana_pid[i_mass]) delete tg_v2_BW_ana_pid[i_mass];
        tg_v2_BW_ana_pid[i_mass] = new TGraph();
        for(Int_t i_pT = 0; i_pT < N_points_BW_ana; i_pT++)
        {
            Double_t pt_BW = i_pT*delta_pT + 0.0;
            if(pt_BW > min_max_pT_range_pid_plot_ana[0][i_mass] && pt_BW < min_max_pT_range_pid_plot_ana[1][i_mass])
            {
                Double_t v2_BW = 0;
                Double_t inv_yield_BW = 0;
                if(id_bw_hypersurface == 1) bw_ana.calc_blastwave_yield_and_v2_fos1(pt_BW, arr_quark_mass_meson[i_mass], T_BW, Rho0_BW, Rho2_BW, RxOverRy_BW, inv_yield_BW, v2_BW);
                if(id_bw_hypersurface == 2) bw_ana.calc_blastwave_yield_and_v2_fos2(pt_BW, arr_quark_mass_meson[i_mass], T_BW, Rho0_BW, Rho2_BW, RxOverRy_BW, inv_yield_BW, v2_BW);
                if(id_bw_hypersurface == 3) bw_ana.calc_blastwave_yield_and_v2_fos3(pt_BW, arr_quark_mass_meson[i_mass], T_BW, Rho0_BW, Rho2_BW, RxOverRy_BW, inv_yield_BW, v2_BW);
                tg_v2_BW_ana_pid[i_mass] ->SetPoint(i_pT,pt_BW,v2_BW);
            }
        }
        //tg_v2_BW_ana_pid[i_mass] -> SetLineColor(arr_color_mass[i_mass]);
        tg_v2_BW_ana_pid[i_mass] -> SetLineColor(arr_color_mass[i_mass]);
        tg_v2_BW_ana_pid[i_mass] -> SetLineWidth(3);
        tg_v2_BW_ana_pid[i_mass] -> SetLineStyle(9);
        //if(!(fCheckBox_pid[i_mass]->GetState() == kButtonDown)) tg_v2_BW_ana_pid[i_mass] -> SetLineStyle(9);
        if(fCheckBox_sel[1]->GetState() == kButtonDown && fCheckBox_pid[i_mass]->GetState() == kButtonDown) tg_v2_BW_ana_pid[i_mass] -> Draw("same L");
    }

    printf("v2 ana plotted \n");
    c_1X1_v2->GetCanvas()->Modified();
    c_1X1_v2->GetCanvas()->Update();
    //------------------------------------------------------



    //------------------------------------------------------
    // dNdpT plots
    c_1X1_dNdpt ->GetCanvas() ->cd();

    for(int i_mass = 0; i_mass < N_masses_all; ++i_mass)
    {
        if(!(fCheckBox_pid_fit_dNdpt[i_mass]->GetState() == kButtonDown)) continue;
        //Double_t delta_pT = (Double_t)((integration_range_pid[i_mass][1] - 0.0)/(Double_t)N_points_BW_ana);
        Double_t delta_pT = (Double_t)((min_max_pT_range_pid_plot_ana[1][i_mass] - min_max_pT_range_pid_plot_ana[0][i_mass])/(Double_t)(N_points_BW_ana-1));
        //Double_t min_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();
        //Double_t max_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();

        //fCanvasB ->GetCanvas()->cd(i_mass + 1);
        if(tg_dNdpT_BW_ana_pid[i_mass]) delete tg_dNdpT_BW_ana_pid[i_mass];
        tg_dNdpT_BW_ana_pid[i_mass] = new TGraph();
        Double_t integral_BW = 0.0;
        for(Int_t i_pT = 0; i_pT < N_points_BW_ana; i_pT++)
        {
            Double_t pt_BW = i_pT*delta_pT + 0.0;
            Double_t v2_BW = 0;
            Double_t inv_yield_BW = 0;
            if(id_bw_hypersurface == 1) bw_ana.calc_blastwave_yield_and_v2_fos1(pt_BW, arr_quark_mass_meson[i_mass], T_BW, Rho0_BW, Rho2_BW, RxOverRy_BW, inv_yield_BW, v2_BW);
            if(id_bw_hypersurface == 2) bw_ana.calc_blastwave_yield_and_v2_fos2(pt_BW, arr_quark_mass_meson[i_mass], T_BW, Rho0_BW, Rho2_BW, RxOverRy_BW, inv_yield_BW, v2_BW);
            if(id_bw_hypersurface == 3) bw_ana.calc_blastwave_yield_and_v2_fos3(pt_BW, arr_quark_mass_meson[i_mass], T_BW, Rho0_BW, Rho2_BW, RxOverRy_BW, inv_yield_BW, v2_BW);
            inv_yield_BW *= pt_BW;
            tg_dNdpT_BW_ana_pid[i_mass] ->SetPoint(i_pT,pt_BW,inv_yield_BW);

            if(pt_BW > integration_range_pid[i_mass][0] && pt_BW < integration_range_pid[i_mass][1])
            {
                integral_BW += inv_yield_BW*delta_pT;
            }
        }
        if(integral_BW > 0.0)
        {
            for(Int_t i_pT = 0; i_pT < tg_dNdpT_BW_ana_pid[i_mass]->GetN(); i_pT++)
            {
                Double_t pt_BW, y_val_BW;
                tg_dNdpT_BW_ana_pid[i_mass] ->GetPoint(i_pT,pt_BW,y_val_BW);
                tg_dNdpT_BW_ana_pid[i_mass] ->SetPoint(i_pT,pt_BW,y_val_BW/integral_BW);
                cout << "i_pT = " << i_pT << ", pt_BW= " << pt_BW << ", y_val_BW= " << y_val_BW << ", y_val_BW/integral_BW= " << y_val_BW/integral_BW << endl;
            }
        }
        cout << "integral_BW" <<integral_BW <<endl;
        tg_dNdpT_BW_ana_pid[i_mass] -> SetLineColor(arr_color_mass[i_mass]);
        tg_dNdpT_BW_ana_pid[i_mass] -> SetLineWidth(3);
        tg_dNdpT_BW_ana_pid[i_mass] -> SetLineStyle(1);
        //if(!(fCheckBox_pid[i_mass]->GetState() == kButtonDown)) tg_dNdpT_BW_ana_pid[i_mass] -> SetLineStyle(9);
        if(fCheckBox_sel[1]->GetState() == kButtonDown && fCheckBox_pid_fit_dNdpt[i_mass]->GetState() == kButtonDown)  tg_dNdpT_BW_ana_pid[i_mass] -> Draw("same L");
    }

    printf("dNdpT ana plotted \n");
    c_1X1_dNdpt->GetCanvas()->Modified();
    c_1X1_dNdpt->GetCanvas()->Update();
    //------------------------------------------------------
}
//______________________________________________________________________________



//______________________________________________________________________________
void TBlastWaveGUI::MakePlotv2()
{
    printf("TBlastWaveGUI::MakePlotv2() \n");

    Int_t id_bw_hypersurface = fCombo ->GetSelected();
    printf("id_bw_hypersurface: %d \n",id_bw_hypersurface);

    Int_t i_Temp   = vec_slider[0]->GetPosition();
    Int_t i_rho_0  = vec_slider[1]->GetPosition();
    Int_t i_rho_a  = vec_slider[2]->GetPosition();
    Int_t i_R_x    = vec_slider[3]->GetPosition();
    Int_t i_fboost = vec_slider[4]->GetPosition();

    //--------------------------------------------------------------------
    Double_t min_max_pT_range_pid_plot_ana[2][N_masses];
    min_max_pT_range_pid_plot_ana[0][0] = 0.1; // pi
    min_max_pT_range_pid_plot_ana[1][0] = 3.5;
    min_max_pT_range_pid_plot_ana[0][1] = 0.1; // K
    min_max_pT_range_pid_plot_ana[1][1] = 4.0;
    min_max_pT_range_pid_plot_ana[0][2] = 0.1; // p
    min_max_pT_range_pid_plot_ana[1][2] = 5.0;
    min_max_pT_range_pid_plot_ana[0][3] = 0.1; // phi
    min_max_pT_range_pid_plot_ana[1][3] = 5.0;
    min_max_pT_range_pid_plot_ana[0][4] = 0.1; // Omega
    min_max_pT_range_pid_plot_ana[1][4] = 6.0;
    min_max_pT_range_pid_plot_ana[0][5] = 0.1; // D0
    min_max_pT_range_pid_plot_ana[1][5] = 12.5;
    min_max_pT_range_pid_plot_ana[0][6] = 0.1; // J/Psi
    min_max_pT_range_pid_plot_ana[1][6] = 15.0;
    min_max_pT_range_pid_plot_ana[0][7] = 0.1; // Upsilon
    min_max_pT_range_pid_plot_ana[1][7] = 15.0;
    min_max_pT_range_pid_plot_ana[0][8] = 0.1; // d
    min_max_pT_range_pid_plot_ana[1][8] = 10.0;
    //--------------------------------------------------------------------

    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_make_plot_v2->ChangeBackground(green);


    //------------------------------------------------------------------------------------------------------------------------------------
    Double_t x_range_plot[3][2] =
    {
        {-0.35,4.3},
        {-0.35,8.2},
        {-0.35,10.7}
    };

    Double_t y_range_plot[3][2] =
    {
        {-0.05,0.34},
        {-0.05,0.34},
        {-0.05,0.34}
    };

    Int_t arr_color_plot[3][4] =
    {
        {kBlack,kBlue+1,kGreen+2,kGreen+2},
        {kBlack,kBlue+1,kGreen+2,kGreen+2},
        {kGreen+2,kRed,kBlue+1,kGreen+2}
    };
    Double_t arr_size_plot[3][4] =
    {
        {1.2,1.3,1.7,1.4},
        {1.2,1.3,1.7,1.4},
        {1.2,1.3,1.5,1.4}
    };
    Int_t arr_marker_styleAB[3][4] =
    {
        {25,46,29,28},
        {25,46,29,28},
        {21,47,30,28}
    };
    Int_t arr_marker_styleAB_x1[3][4] =
    {
        {21,47,29,34},
        {21,47,29,34},
        {21,47,29,34}
    };

    // pi, K, p, phi, Omega, D0, J/Psi, Upsilon, d
    Int_t arr_plot_orderAB[3][3] = // [pad][particle]
    {
        {0,8,3}, // pi, d, phi
        {2,4,6}, // p, Omega, J/Psi
        {5,7,1}  // D0, Upsilon, K
    };

    if(!c_3X1)
    {
        c_3X1 = new TCanvas("c_3X1","c_3X1",100,200,1400,600);
        c_3X1->SetFillColor(10);
        c_3X1->SetTopMargin(0.05);
        c_3X1->SetBottomMargin(0.22);
        c_3X1->SetRightMargin(0.05);
        c_3X1->SetLeftMargin(0.22);
        c_3X1->SetLogy(0);
        c_3X1->Divide(3,1,0.0,0.0,10);
    }


    for(Int_t iPad = 0; iPad < 3; iPad++)
    {
        Double_t scaling_factor_3X1 = 1.0;
        Double_t Label_size_3X1     = 0.07;
        if(iPad == 0)
        {
            scaling_factor_3X1 = 0.78;
        }
        c_3X1->cd(iPad+1)->SetTicks(1,1);
        c_3X1->cd(iPad+1)->SetGrid(0,0);
        c_3X1->cd(iPad+1)->SetFillColor(10);
        c_3X1->cd(iPad+1)->SetRightMargin(0.01);
        c_3X1->cd(iPad+1)->SetTopMargin(0.01);

        HistName = "h_frame_3X1_";
        HistName += iPad;

        if(h_frame_3X1[iPad]) delete h_frame_3X1[iPad];
        h_frame_3X1[iPad] = c_3X1->cd(iPad+1)->DrawFrame(x_range_plot[iPad][0],y_range_plot[iPad][0],x_range_plot[iPad][1],y_range_plot[iPad][1],HistName.Data());
        h_frame_3X1[iPad]->SetStats(0);
        h_frame_3X1[iPad]->SetTitle("");
        h_frame_3X1[iPad]->GetXaxis()->SetTitleOffset(0.95/scaling_factor_3X1);
        h_frame_3X1[iPad]->GetYaxis()->SetTitleOffset(1.6);
        h_frame_3X1[iPad]->GetXaxis()->SetLabelOffset(0.0);
        if(iPad == 0) h_frame_3X1[iPad]->GetXaxis()->SetLabelOffset(0.014);
        //if(iPad == 0) h_frame_3X1[iPad]->GetXaxis()->SetLabelOffset(0.0*scaling_factor_3X1);
        h_frame_3X1[iPad]->GetYaxis()->SetLabelOffset(0.01*scaling_factor_3X1);
        h_frame_3X1[iPad]->GetXaxis()->SetLabelSize(Label_size_3X1*scaling_factor_3X1);
        h_frame_3X1[iPad]->GetYaxis()->SetLabelSize(Label_size_3X1*scaling_factor_3X1);
        h_frame_3X1[iPad]->GetXaxis()->SetTitleSize(Label_size_3X1*scaling_factor_3X1);
        h_frame_3X1[iPad]->GetYaxis()->SetTitleSize(Label_size_3X1*scaling_factor_3X1);
        h_frame_3X1[iPad]->GetXaxis()->SetNdivisions(505,'N');
        h_frame_3X1[iPad]->GetYaxis()->SetNdivisions(505,'N');
        h_frame_3X1[iPad]->GetXaxis()->CenterTitle();
        h_frame_3X1[iPad]->GetYaxis()->CenterTitle();
        if(iPad == 1) h_frame_3X1[iPad]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        else h_frame_3X1[iPad]->GetXaxis()->SetTitle("");
        h_frame_3X1[iPad]->GetYaxis()->SetTitle("v_{2}");
        h_frame_3X1[iPad]->GetXaxis()->SetRangeUser(x_range_plot[iPad][0],x_range_plot[iPad][1]);
        h_frame_3X1[iPad]->GetYaxis()->SetRangeUser(y_range_plot[iPad][0],y_range_plot[iPad][1]);

        if(TL_line_base_plot[iPad]) delete TL_line_base_plot[iPad];
        TL_line_base_plot[iPad] = PlotLine(x_range_plot[iPad][0],x_range_plot[iPad][1],0.0,0.0,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)

        const Int_t N_points_BW_ana = 60;
        printf("T_BW_fit_ana: %4.3f \n",T_BW_fit_ana);

        Int_t i_particle_use = 0;
        for(Int_t i_mass_loop = 0; i_mass_loop < 3; i_mass_loop++)
        {
            Int_t i_mass = arr_plot_orderAB[iPad][i_mass_loop];
            if(fCheckBox_pid_plot[i_mass]->GetState() == kButtonDown)
            {
                //printf("Draw i_mass: %d \n", i_mass);
                tgae_v2_vs_pT_mesons_data_copy[i_mass] ->SetMarkerColor(kGray+1);
                tgae_v2_vs_pT_mesons_data_copy[i_mass] ->SetLineColor(kGray);
                tgae_v2_vs_pT_mesons_data_copy[i_mass] ->SetMarkerSize(arr_size_plot[iPad][i_mass_loop]*1.4);
                tgae_v2_vs_pT_mesons_data_copy[i_mass] ->SetMarkerStyle(arr_marker_styleAB_x1[iPad][i_mass_loop]);
                tgae_v2_vs_pT_mesons_data_copy[i_mass] ->Draw("same P");

                tgae_v2_vs_pT_mesons_data_copyB[i_mass] ->SetMarkerColor(kWhite);
                tgae_v2_vs_pT_mesons_data_copyB[i_mass] ->SetLineColor(kGray);
                tgae_v2_vs_pT_mesons_data_copyB[i_mass] ->SetMarkerSize(arr_size_plot[iPad][i_mass_loop]*1.0);
                tgae_v2_vs_pT_mesons_data_copyB[i_mass] ->SetMarkerStyle(arr_marker_styleAB_x1[iPad][i_mass_loop]);
                tgae_v2_vs_pT_mesons_data_copyB[i_mass] ->Draw("same P");

                tgae_v2_vs_pT_mesons_data[i_mass] ->SetMarkerColor(arr_color_plot[iPad][i_mass_loop]);
                tgae_v2_vs_pT_mesons_data[i_mass] ->SetLineColor(kGray);
                tgae_v2_vs_pT_mesons_data[i_mass] ->SetMarkerSize(arr_size_plot[iPad][i_mass_loop]);
                tgae_v2_vs_pT_mesons_data[i_mass] ->SetMarkerStyle(arr_marker_styleAB[iPad][i_mass_loop]);
                tgae_v2_vs_pT_mesons_data[i_mass] ->Draw("same P");
                // OK


                tg_label_plot[0][i_mass] ->SetMarkerColor(kGray+1);
                tg_label_plot[0][i_mass] ->SetMarkerSize(arr_size_plot[iPad][i_mass_loop]*1.4);
                tg_label_plot[0][i_mass] ->SetMarkerStyle(arr_marker_styleAB_x1[iPad][i_mass_loop]);

                tg_label_plot[1][i_mass] ->SetMarkerColor(kWhite);
                tg_label_plot[1][i_mass] ->SetMarkerSize(arr_size_plot[iPad][i_mass_loop]*1.0);
                tg_label_plot[1][i_mass] ->SetMarkerStyle(arr_marker_styleAB_x1[iPad][i_mass_loop]);

                tg_label_plot[2][i_mass] ->SetMarkerColor(arr_color_plot[iPad][i_mass_loop]);
                tg_label_plot[2][i_mass] ->SetMarkerSize(arr_size_plot[iPad][i_mass_loop]);
                tg_label_plot[2][i_mass] ->SetMarkerStyle(arr_marker_styleAB[iPad][i_mass_loop]);

                Double_t offset_x_pad[3] = {0.36,0.84,0.84};
                Double_t offset_y_pad[3] = {0.91,0.91,0.91};

                Double_t x_pos_legend = offset_x_pad[iPad];
                Double_t y_pos_legend = offset_y_pad[iPad]  - i_particle_use*0.055;

                Double_t x_user, y_user;
                Double_t x_NDC = offset_x_pad[iPad];
                Double_t y_NDC = offset_y_pad[iPad]  - i_particle_use*0.055;
                get_user_from_NDC((TPad*)c_3X1->cd(iPad+1),x_NDC-0.026,y_NDC+0.013,x_user,y_user);
                printf("iPad: %d, (x,y): {%4.3f, %4.3f} \n",iPad,x_user,y_user);

                tg_label_plot[0][i_mass] ->SetPoint(0,x_user,y_user);
                tg_label_plot[1][i_mass] ->SetPoint(0,x_user,y_user);
                tg_label_plot[2][i_mass] ->SetPoint(0,x_user,y_user);
                tg_label_plot[0][i_mass] ->Draw("same P");
                tg_label_plot[1][i_mass] ->Draw("same P");
                tg_label_plot[2][i_mass] ->Draw("same P");


                if(TLatex_legend_v2_plot[i_mass]) delete TLatex_legend_v2_plot[i_mass];
                TLatex_legend_v2_plot[i_mass] = plotTopLegend((char*)label_pid_spectra[i_mass].Data(),x_pos_legend,y_pos_legend,0.065*scaling_factor_3X1,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
                i_particle_use++;


                //-----------------------------
                Double_t min_val_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();
                Double_t max_val_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();

                Double_t delta_pT = (Double_t)((min_max_pT_range_pid_plot_ana[1][i_mass] - min_max_pT_range_pid_plot_ana[0][i_mass])/(Double_t)(N_points_BW_ana-1));
                //Double_t min_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();
                //Double_t max_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();


                if(tg_v2_BW_ana_pid_plot[i_mass]) delete tg_v2_BW_ana_pid_plot[i_mass];
                if(tg_v2_BW_ana_pid_plot_range[i_mass]) delete tg_v2_BW_ana_pid_plot_range[i_mass];
                tg_v2_BW_ana_pid_plot[i_mass]       = new TGraph();
                tg_v2_BW_ana_pid_plot_range[i_mass] = new TGraph();

                Int_t point_counter       = 0;
                Int_t point_counter_range = 0;
                for(Int_t i_pT = 0; i_pT < N_points_BW_ana; i_pT++)
                {
                    Double_t pt_BW = i_pT*delta_pT + 0.0;
                    if(pt_BW > min_max_pT_range_pid_plot_ana[0][i_mass] && pt_BW < min_max_pT_range_pid_plot_ana[1][i_mass])
                    {
                        Double_t v2_BW = 0;
                        Double_t inv_yield_BW = 0;

                        if(id_bw_hypersurface == 1) bw_ana.calc_blastwave_yield_and_v2_fos1(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
                        if(id_bw_hypersurface == 2) bw_ana.calc_blastwave_yield_and_v2_fos2(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
                        if(id_bw_hypersurface == 3) bw_ana.calc_blastwave_yield_and_v2_fos3(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
                        tg_v2_BW_ana_pid_plot[i_mass]       ->SetPoint(point_counter,pt_BW,v2_BW);
                        point_counter++;
                        if(pt_BW > min_val_pT && pt_BW < max_val_pT)
                        {
                            tg_v2_BW_ana_pid_plot_range[i_mass] ->SetPoint(point_counter_range,pt_BW,v2_BW);
                            point_counter_range++;
                        }
                    }
                }
                tg_v2_BW_ana_pid_plot[i_mass] -> SetLineColorAlpha(kGray+1,0.5);
                tg_v2_BW_ana_pid_plot[i_mass] -> SetLineWidth(2);
                tg_v2_BW_ana_pid_plot[i_mass] -> SetLineStyle(9);

                tg_v2_BW_ana_pid_plot_range[i_mass] -> SetLineColor(arr_color_plot[iPad][i_mass_loop]);
                tg_v2_BW_ana_pid_plot_range[i_mass] -> SetLineWidth(3);
                if(fCheckBox_pid[i_mass]->GetState() == kButtonDown) tg_v2_BW_ana_pid_plot_range[i_mass] -> SetLineStyle(1);
                else tg_v2_BW_ana_pid_plot_range[i_mass] -> SetLineStyle(9);

                Double_t chi2 = 0.0;
                if((fCheckBox_pid[i_mass]->GetState() == kButtonDown))
                {
                    //------------------------------------------------------------
                    // Calculate chi2
                    Int_t N_points_data = tgae_v2_vs_pT_mesons_data[i_mass]->GetN();
                    Int_t N_points_data_used = 0;
                    for(int i_point = 0; i_point < N_points_data; ++i_point)
                    {
                        double v2_data = 0.0;
                        double pt_data    = 0.0;

                        tgae_v2_vs_pT_mesons_data[i_mass]->GetPoint(i_point,pt_data,v2_data);
                        if(pt_data < min_val_pT || pt_data > max_val_pT) continue;
                        double v2_err = tgae_v2_vs_pT_mesons_data[i_mass]->GetErrorYhigh(i_point);

                        Double_t y_val_BW = tg_v2_BW_ana_pid_plot[i_mass]->Eval(pt_data);
                        chi2 += TMath::Power((v2_data - y_val_BW)/v2_err,2);
                        N_points_data_used++;

                    }
                    if(N_points_data > 0) chi2 /= (Double_t)N_points_data_used;
                    //------------------------------------------------------------
                }
                //-----------------------------

                Double_t x_offset_chi2 = 0.3;
                Int_t align_chi2 = 32;
                if(iPad > 0)
                {
                    x_offset_chi2 = -0.3;
                    align_chi2    = 1;
                }
                if(fCheckBox_sel[2]->GetState() == kButtonDown)
                {
                    HistName = "#chi^{2} = ";
                    //sprintf(NoP,"%4.2f",best_individual_chi2_per_point_ana[i_mass][0]);
                    sprintf(NoP,"%4.2f",chi2);
                    HistName += NoP;
                    plotTopLegend((char*)HistName.Data(),x_pos_legend+x_offset_chi2,y_pos_legend,0.065*scaling_factor_3X1,kBlack,0.0,42,1,align_chi2); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
                }


                if(fCheckBox_sel[0]->GetState() == kButtonDown)
                {
                    tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineWidth(3);
                    tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineColor(arr_color_plot[iPad][i_mass_loop]);
                    tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->Draw("same L hist");
                }

            }

        }


        if(T_BW_fit_ana > 0.0)
        {
            for(int i_mass_loop = 0; i_mass_loop < 3; ++i_mass_loop)
            {
                Int_t i_mass = arr_plot_orderAB[iPad][i_mass_loop];
                if((fCheckBox_pid_plot[i_mass]->GetState() == kButtonDown))
                {
                    if(fCheckBox_sel[1]->GetState() == kButtonDown)
                    {
                        tg_v2_BW_ana_pid_plot[i_mass]       -> Draw("same L");
                        tg_v2_BW_ana_pid_plot_range[i_mass] -> Draw("same L");
                    }
                }
            }
        }

        /*
        if(T_BW_fit_ana > 0.0)
        {
            for(int i_mass_loop = 0; i_mass_loop < 3; ++i_mass_loop)
            {
                Int_t i_mass = arr_plot_orderAB[iPad][i_mass_loop];

                Double_t min_val_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();
                Double_t max_val_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();

                Double_t delta_pT = (Double_t)((min_max_pT_range_pid_plot_ana[1][i_mass] - min_max_pT_range_pid_plot_ana[0][i_mass])/(Double_t)N_points_BW_ana);
                //Double_t min_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();
                //Double_t max_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();


                if(tg_v2_BW_ana_pid_plot[i_mass]) delete tg_v2_BW_ana_pid_plot[i_mass];
                if(tg_v2_BW_ana_pid_plot_range[i_mass]) delete tg_v2_BW_ana_pid_plot_range[i_mass];
                tg_v2_BW_ana_pid_plot[i_mass]       = new TGraph();
                tg_v2_BW_ana_pid_plot_range[i_mass] = new TGraph();
                for(Int_t i_pT = 0; i_pT < N_points_BW_ana; i_pT++)
                {
                    Double_t pt_BW = i_pT*delta_pT + 0.0;
                    if(pt_BW > min_max_pT_range_pid_plot_ana[0][i_mass] && pt_BW < min_max_pT_range_pid_plot_ana[1][i_mass])
                    {
                        Double_t v2_BW = 0;
                        Double_t inv_yield_BW = 0;

                        if(id_bw_hypersurface == 1) bw_ana.calc_blastwave_yield_and_v2_fos1(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
                        if(id_bw_hypersurface == 2) bw_ana.calc_blastwave_yield_and_v2_fos2(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
                        if(id_bw_hypersurface == 3) bw_ana.calc_blastwave_yield_and_v2_fos3(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
                        tg_v2_BW_ana_pid_plot[i_mass]       ->SetPoint(i_pT,pt_BW,v2_BW);
                        if(pt_BW > min_val_pT && pt_BW < max_val_pT)
                        {
                            tg_v2_BW_ana_pid_plot_range[i_mass] ->SetPoint(i_pT,pt_BW,v2_BW);
                        }
                    }
                }
                tg_v2_BW_ana_pid_plot[i_mass] -> SetLineColorAlpha(kGray+1,0.5);
                tg_v2_BW_ana_pid_plot[i_mass] -> SetLineWidth(2);
                tg_v2_BW_ana_pid_plot[i_mass] -> SetLineStyle(9);

                tg_v2_BW_ana_pid_plot_range[i_mass] -> SetLineColor(arr_color_plot[iPad][i_mass_loop]);
                tg_v2_BW_ana_pid_plot_range[i_mass] -> SetLineWidth(3);
                tg_v2_BW_ana_pid_plot_range[i_mass] -> SetLineStyle(9);
                if((fCheckBox_pid[i_mass]->GetState() == kButtonDown))
                {

                    //------------------------------------------------------------
                    // Calculate chi2
                    Int_t N_points_data = tgae_v2_vs_pT_mesons_data[i_mass]->GetN();
                    for(int i_point = 0; i_point < N_points_data; ++i_point)
                    {
                        double v2_data = 0.0;
                        double pt_data    = 0.0;

                        tgae_v2_vs_pT_mesons_data[i_mass]->GetPoint(i_point,pt_data,v2_data);
                        double v2_err = tgae_v2_vs_pT_mesons_data[i_mass]->GetErrorYhigh(i_point);

                        Double_t y_val_BW = [i_mass] tg_v2_BW_ana_pid_plot->Eval(pt_data);
                        chi2 += TMath::Power((v2_data - y_val_BW)/v2_err,2);

                    }
                    if(N_points_data > 0) chi2 /= (Double_t)N_points_data;
                    //------------------------------------------------------------

                    if(fCheckBox_sel[1]->GetState() == kButtonDown)
                    {
                        tg_v2_BW_ana_pid_plot[i_mass]       -> Draw("same L");
                        tg_v2_BW_ana_pid_plot_range[i_mass] -> Draw("same L");
                    }
                }
            }
        }
        */

        if(iPad == 0)
        {
            tg_leg ->SetLineColor(kBlack);
            tg_leg ->SetLineWidth(3);
            tg_leg ->SetLineStyle(1);

            leg_v2_vs_pT_C = new TLegend(0.68,0.4,0.98,0.5); // x1,y1,x2,y2
            leg_v2_vs_pT_C->SetBorderSize(0);
            leg_v2_vs_pT_C->SetFillColor(0);
            leg_v2_vs_pT_C->SetTextSize(Label_size_3X1*scaling_factor_3X1);
            leg_v2_vs_pT_C->AddEntry((TGraph*)tg_leg->Clone(),"fit","l");
            tg_leg ->SetLineStyle(9);
            leg_v2_vs_pT_C->AddEntry((TGraph*)tg_leg->Clone(),"prediction","l");
            leg_v2_vs_pT_C->Draw();
        }

    } // end of iPad loop

    c_3X1 ->Modified();
    c_3X1 ->Update();
    printf("v2 plotted \n");
    //------------------------------------------------------------------------------------------------------------------------------------

}
//______________________________________________________________________________



//______________________________________________________________________________
void TBlastWaveGUI::MakePlotdNdpT()
{
    printf("TBlastWaveGUI::MakePlotdNdpT() \n");

    Int_t id_bw_hypersurface = fCombo ->GetSelected();
    printf("id_bw_hypersurface: %d \n",id_bw_hypersurface);

    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_make_plot_dNdpT->ChangeBackground(green);

    Int_t i_Temp   = vec_slider[0]->GetPosition();
    Int_t i_rho_0  = vec_slider[1]->GetPosition();
    Int_t i_rho_a  = vec_slider[2]->GetPosition();
    Int_t i_R_x    = vec_slider[3]->GetPosition();
    Int_t i_fboost = vec_slider[4]->GetPosition();

    //--------------------------------------------------------------------
    Double_t min_max_pT_range_pid_plot_ana[2][N_masses];
    min_max_pT_range_pid_plot_ana[0][0] = 0.1; // pi
    min_max_pT_range_pid_plot_ana[1][0] = 4.5;
    min_max_pT_range_pid_plot_ana[0][1] = 0.1; // K
    min_max_pT_range_pid_plot_ana[1][1] = 4.4;
    min_max_pT_range_pid_plot_ana[0][2] = 0.1; // p
    min_max_pT_range_pid_plot_ana[1][2] = 4.5;
    min_max_pT_range_pid_plot_ana[0][3] = 0.1; // phi
    min_max_pT_range_pid_plot_ana[1][3] = 4.5;
    min_max_pT_range_pid_plot_ana[0][4] = 0.1; // Omega
    min_max_pT_range_pid_plot_ana[1][4] = 13.5;
    min_max_pT_range_pid_plot_ana[0][5] = 0.1; // D0
    min_max_pT_range_pid_plot_ana[1][5] = 13.5;
    min_max_pT_range_pid_plot_ana[0][6] = 0.1; // J/Psi
    min_max_pT_range_pid_plot_ana[1][6] = 13.5;
    min_max_pT_range_pid_plot_ana[0][7] = 0.1; // Upsilon
    min_max_pT_range_pid_plot_ana[1][7] = 13.5;
    min_max_pT_range_pid_plot_ana[0][8] = 0.1; // Upsilon
    min_max_pT_range_pid_plot_ana[1][8] = 13.5;
    //--------------------------------------------------------------------

    if(!c_3X3)
    {
        c_3X3 = new TCanvas("c_3X3","c_3X3",500,20,950,900);
        c_3X3->SetTopMargin(0.15);
        c_3X3->SetBottomMargin(0.25);
        c_3X3->SetRightMargin(0.15);
        c_3X3->SetLeftMargin(0.25);

        c_3X3->SetLogy(0);
        c_3X3->Divide(3,3,0.0,0.0); // x divide, y divide, x margin, y margin
      }

    const Int_t N_points_BW_ana = 60; // 55
    // 0   1   2  3    4      5    6      7      8
    // pi, K, p, phi, Omega, D0, J/Psi, Upsilon, d
    Int_t arr_plot_order[9] = {0,3,8,1,4,6,2,5,7};
    Double_t x_axis_range[9] = {3.5,5.5,10.5,3.5,5.5,10.5,3.5,5.5,17.5};

    Double_t y_max_plot = 2.2;
    Double_t Label_size_3X3     = 0.1;

    for(Int_t iPad = 0; iPad < N_masses; iPad++)
    {
        Int_t i_mass = arr_plot_order[iPad];

        if(h_frame_3X3[iPad]) delete h_frame_3X3[iPad];

        Double_t scaling_factor_3X3 = 1.0;
        Double_t offset_3X3         = 0.0;
        Double_t Label_size_2X4     = 0.1;
        if(iPad == 6)
        {
            scaling_factor_3X3 = 0.77;
            offset_3X3         = 0.015;
        }
        if(iPad > 6)
        {
            //scaling_factor_3X3 = 0.88;
            scaling_factor_3X3 = 0.94;
            offset_3X3         = 0.005;
        }


        c_3X3->cd(iPad+1)->SetTicks(1,1);
        c_3X3->cd(iPad+1)->SetGrid(0,0);
        c_3X3->cd(iPad+1)->SetFillColor(10);
        c_3X3->cd(iPad+1)->SetRightMargin(0.01);
        c_3X3->cd(iPad+1)->SetTopMargin(0.01);
        HistName = "h_frame_3X3_";
        HistName += iPad;
        //h_frame_3X3[iPad] = c_3X3->cd(iPad+1)->DrawFrame(-0.1,-0.1,min_max_pT_range_pid_plot_ana[1][i_mass],y_max_plot,HistName.Data());
        h_frame_3X3[iPad] = c_3X3->cd(iPad+1)->DrawFrame(-0.5,-0.1,x_axis_range[iPad],y_max_plot,HistName.Data());
        h_frame_3X3[iPad]->SetStats(0);
        h_frame_3X3[iPad]->SetTitle("");
        h_frame_3X3[iPad]->SetTickLength(0.04,"X");
        h_frame_3X3[iPad]->SetTickLength(0.04,"Y");
        h_frame_3X3[iPad]->GetXaxis()->SetTitleOffset(0.9/scaling_factor_3X3);
        h_frame_3X3[iPad]->GetYaxis()->SetTitleOffset(1.2);
        h_frame_3X3[iPad]->GetXaxis()->SetLabelOffset(0.0*scaling_factor_3X3 + offset_3X3);
        if(iPad == 0) h_frame_3X3[iPad]->GetXaxis()->SetLabelOffset(0.0*scaling_factor_3X3);
        h_frame_3X3[iPad]->GetYaxis()->SetLabelOffset(0.01*scaling_factor_3X3);
        h_frame_3X3[iPad]->GetXaxis()->SetLabelSize(Label_size_3X3*scaling_factor_3X3);
        h_frame_3X3[iPad]->GetYaxis()->SetLabelSize(Label_size_3X3*scaling_factor_3X3);
        h_frame_3X3[iPad]->GetXaxis()->SetTitleSize(Label_size_3X3*scaling_factor_3X3);
        h_frame_3X3[iPad]->GetYaxis()->SetTitleSize(Label_size_3X3*scaling_factor_3X3);
        h_frame_3X3[iPad]->GetXaxis()->SetNdivisions(505,'N');
        h_frame_3X3[iPad]->GetYaxis()->SetNdivisions(505,'N');
        h_frame_3X3[iPad]->GetXaxis()->CenterTitle();
        h_frame_3X3[iPad]->GetYaxis()->CenterTitle();
        if(iPad == 7) h_frame_3X3[iPad]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        else h_frame_3X3[iPad]->GetXaxis()->SetTitle("");
        if(iPad == 3) h_frame_3X3[iPad]->GetYaxis()->SetTitle("dN/dp_{T} (GeV/c)^{-1}");
        else h_frame_3X3[iPad]->GetYaxis()->SetTitle("");
        if(iPad < 3) h_frame_3X3[iPad]->GetYaxis()->SetRangeUser(-0.1,2.2);
        if(iPad >= 3 && iPad < 6) h_frame_3X3[iPad]->GetYaxis()->SetRangeUser(-0.05,1.2);
        if(iPad >= 6) h_frame_3X3[iPad]->GetYaxis()->SetRangeUser(-0.04,0.85);

        if(tg_dNdpT_BW_ana_pid_plot[i_mass]) delete tg_dNdpT_BW_ana_pid_plot[i_mass];
        tg_dNdpT_BW_ana_pid_plot[i_mass] = new TGraph();
        if(tg_dNdpT_BW_ana_pid_plot_range[i_mass]) delete tg_dNdpT_BW_ana_pid_plot_range[i_mass];
        tg_dNdpT_BW_ana_pid_plot_range[i_mass] = new TGraph();
        Double_t integral_BW = 0.0;

        Double_t min_val_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();
        Double_t max_val_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();

        tgae_dN_dpT_mesons_data_A[i_mass] ->SetMarkerColor(kBlack);
        tgae_dN_dpT_mesons_data_A[i_mass] ->SetMarkerStyle(20);
        tgae_dN_dpT_mesons_data_A[i_mass] ->SetLineColor(kGray);
        tgae_dN_dpT_mesons_data_A[i_mass] ->SetMarkerSize(1.15);

        tgae_dN_dpT_mesons_data_B[i_mass] ->SetMarkerColor(kWhite);
        tgae_dN_dpT_mesons_data_B[i_mass] ->SetMarkerStyle(20);
        tgae_dN_dpT_mesons_data_B[i_mass] ->SetLineColor(kGray);
        tgae_dN_dpT_mesons_data_B[i_mass] ->SetMarkerSize(0.9);

        tgae_dN_dpT_mesons_data[i_mass] ->SetMarkerColor(kGray);
        tgae_dN_dpT_mesons_data[i_mass] ->SetMarkerStyle(20);
        tgae_dN_dpT_mesons_data[i_mass] ->SetLineColor(kGray);
        tgae_dN_dpT_mesons_data[i_mass] ->SetMarkerSize(0.9);

        if(tgae_dN_dpT_mesons_data[i_mass] && i_mass != 6)
        {
            tgae_dN_dpT_mesons_data_A[i_mass] ->Draw("same P");
            tgae_dN_dpT_mesons_data_B[i_mass] ->Draw("same P");
            tgae_dN_dpT_mesons_data[i_mass]   ->Draw("same P");
        }

        if(fCheckBox_sel[0]->GetState() == kButtonDown
           && h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]
          )
        {
            h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineColor(kBlue); // blast wave MC
            h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineWidth(5); // blast wave MC
            h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->DrawCopy("same hist L"); // blast wave MC
        }


        Double_t chi2 = 0.0;
        if(T_BW_fit_ana > 0.0)
        {
            //--------------------------------
            Double_t scale_factor_ana  = 0.0;

            vector< vector<Double_t> > vec_data_BW;
            vec_data_BW.resize(4); // data, BW, err, pT

            for(int i_point = 0; i_point < tgae_dN_dpT_mesons_data[i_mass]->GetN(); ++i_point)
            {
                double dNdpT_data = 0.0;
                double pt_data    = 0.0;

                tgae_dN_dpT_mesons_data[i_mass]->GetPoint(i_point,pt_data,dNdpT_data);
                if(dNdpT_data <= 0.0) continue;
                double dNdpT_err = tgae_dN_dpT_mesons_data[i_mass]->GetErrorYhigh(i_point);

                Double_t x_err_low_data  = fabs(tgae_dN_dpT_mesons_data[i_mass] ->GetErrorXlow(i_point));
                Double_t x_err_high_data = fabs(tgae_dN_dpT_mesons_data[i_mass] ->GetErrorXhigh(i_point));
                Double_t bin_width = x_err_low_data + x_err_high_data;

                double v2_BW = 0;
                double inv_yield_BW = 0;

                // blast wave parameters
                const double pt_BW = pt_data;         // in GeV

                if(id_bw_hypersurface == 1) bw_ana.calc_blastwave_yield_and_v2_fos1(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
                if(id_bw_hypersurface == 2) bw_ana.calc_blastwave_yield_and_v2_fos2(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
                if(id_bw_hypersurface == 3) bw_ana.calc_blastwave_yield_and_v2_fos3(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
                inv_yield_BW *= pt_BW;

                vec_data_BW[0].push_back(dNdpT_data);
                vec_data_BW[1].push_back(inv_yield_BW);
                vec_data_BW[2].push_back(dNdpT_err);
                vec_data_BW[3].push_back(pt_BW);
            }

            Double_t min_val_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();
            Double_t max_val_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();

            scale_factor_ana = get_norm_scaling_factor_calc(vec_data_BW,min_val_pT,max_val_pT);
            //--------------------------------




            Double_t delta_pT = (Double_t)((min_max_pT_range_pid_plot_ana[1][i_mass] - min_max_pT_range_pid_plot_ana[0][i_mass])/(Double_t)(N_points_BW_ana-1));

            Int_t i_pT_range = 0;
            for(Int_t i_pT = 0; i_pT < N_points_BW_ana; i_pT++)
            {
                Double_t pt_BW = i_pT*delta_pT + 0.0;
                Double_t v2_BW = 0;
                Double_t inv_yield_BW = 0;
                if(id_bw_hypersurface == 1) bw_ana.calc_blastwave_yield_and_v2_fos1(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
                if(id_bw_hypersurface == 2) bw_ana.calc_blastwave_yield_and_v2_fos2(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
                if(id_bw_hypersurface == 3) bw_ana.calc_blastwave_yield_and_v2_fos3(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
                inv_yield_BW *= pt_BW;
                tg_dNdpT_BW_ana_pid_plot[i_mass] ->SetPoint(i_pT,pt_BW,inv_yield_BW);

                if(pt_BW >= min_val_pT && pt_BW <= max_val_pT)
                {
                    tg_dNdpT_BW_ana_pid_plot_range[i_mass] ->SetPoint(i_pT_range,pt_BW,inv_yield_BW);
                    i_pT_range++;
                }

                if(pt_BW > integration_range_pid[i_mass][0] && pt_BW < integration_range_pid[i_mass][1])
                {
                    integral_BW += inv_yield_BW*delta_pT;
                    //if(i_mass == 4) printf("pt_BW: %4.3f, inv_yield_BW: %4.12f, delta_pT: %4.3f, integral_BW: %4.12f \n",pt_BW,inv_yield_BW,delta_pT,integral_BW);
                }
            }
            if(integral_BW > 0.0)
            {
                for(Int_t i_pT = 0; i_pT < tg_dNdpT_BW_ana_pid_plot[i_mass]->GetN(); i_pT++)
                {
                    Double_t pt_BW, y_val_BW;
                    tg_dNdpT_BW_ana_pid_plot[i_mass] ->GetPoint(i_pT,pt_BW,y_val_BW);
                    //tg_dNdpT_BW_ana_pid_plot[i_mass] ->SetPoint(i_pT,pt_BW,y_val_BW/integral_BW);
                    tg_dNdpT_BW_ana_pid_plot[i_mass] ->SetPoint(i_pT,pt_BW,y_val_BW*scale_factor_ana);
                }
                for(Int_t i_pT = 0; i_pT < tg_dNdpT_BW_ana_pid_plot_range[i_mass]->GetN(); i_pT++)
                {
                    Double_t pt_BW, y_val_BW;
                    tg_dNdpT_BW_ana_pid_plot_range[i_mass] ->GetPoint(i_pT,pt_BW,y_val_BW);
                    //tg_dNdpT_BW_ana_pid_plot_range[i_mass] ->SetPoint(i_pT,pt_BW,y_val_BW/integral_BW);
                    tg_dNdpT_BW_ana_pid_plot_range[i_mass] ->SetPoint(i_pT,pt_BW,y_val_BW*scale_factor_ana);
                }
            }



            tg_dNdpT_BW_ana_pid_plot[i_mass] ->SetLineColor(kGray+2);
            tg_dNdpT_BW_ana_pid_plot[i_mass] ->SetLineWidth(4);
            tg_dNdpT_BW_ana_pid_plot[i_mass] ->SetLineStyle(9);

            tg_dNdpT_BW_ana_pid_plot_range[i_mass] ->SetLineColor(kRed);
            tg_dNdpT_BW_ana_pid_plot_range[i_mass] ->SetLineWidth(5);
            tg_dNdpT_BW_ana_pid_plot_range[i_mass] ->SetLineStyle(1);

            //if(!(fCheckBox_pid[i_mass]->GetState() == kButtonDown)) tg_dNdpT_BW_ana_pid_plot[i_mass] -> SetLineStyle(9);
            if(fCheckBox_sel[1]->GetState() == kButtonDown)
            {
                tg_dNdpT_BW_ana_pid_plot[i_mass]       ->DrawClone("same L");

                //------------------------------------------------------------
                // Calculate chi2
                Int_t N_points_data = tgae_dN_dpT_mesons_data[i_mass]->GetN();
                Int_t N_points_data_used = 0;
                for(int i_point = 0; i_point < N_points_data; ++i_point)
                {
                    double dNdpT_data = 0.0;
                    double pt_data    = 0.0;

                    tgae_dN_dpT_mesons_data[i_mass]->GetPoint(i_point,pt_data,dNdpT_data);
                    //printf("pt_data: %4.3f, min_val_pT: %4.3f, max_val_pT: %4.3f \n",pt_data,min_val_pT,max_val_pT);
                    if(pt_data < min_val_pT || pt_data > max_val_pT) continue;
                    double dNdpT_err = tgae_dN_dpT_mesons_data[i_mass]->GetErrorYhigh(i_point);

                    Double_t y_val_BW = tg_dNdpT_BW_ana_pid_plot[i_mass] ->Eval(pt_data);
                    chi2 += TMath::Power((dNdpT_data - y_val_BW)/dNdpT_err,2);
                    N_points_data_used++;

                    //printf("iPad: %d, i_mass: %d, i_point: %d, y_data: %4.3f, y_BW: %4.3f, err: %4.3f, chi2: %4.3f \n",iPad,i_mass,i_point,dNdpT_data,y_val_BW,dNdpT_err,chi2);

                }
                if(N_points_data > 0) chi2 /= (Double_t)N_points_data_used;
                //------------------------------------------------------------

                if((fCheckBox_pid[i_mass]->GetState() == kButtonDown))
                {
                    tg_dNdpT_BW_ana_pid_plot_range[i_mass] ->DrawClone("same L");
                }
            }
        }

        Double_t x_pos_offset = -0.05;
        if(iPad%3 == 0) x_pos_offset = -0.01;
        Double_t x_pos_legend = 0.92 + x_pos_offset;
        Double_t y_pos_legend = 0.88;
        if(TLatex_legend_dNdpT_plot[i_mass]) delete TLatex_legend_dNdpT_plot[i_mass];
        TLatex_legend_dNdpT_plot[i_mass] = plotTopLegend((char*)label_pid_spectra[i_mass].Data(),x_pos_legend,y_pos_legend,0.1*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        if(fCheckBox_sel[2]->GetState() == kButtonDown)
        {
            HistName = "#chi^{2} = ";
            //sprintf(NoP,"%4.2f",best_individual_chi2_per_point_ana[i_mass][1]);
            sprintf(NoP,"%4.2f",chi2);
            HistName += NoP;
        }

        Double_t legend_text_size = 0.08;
        if(iPad == 0) // pi
        {
            plotTopLegend((char*)"30-40%",x_pos_legend,y_pos_legend-0.08*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"2.76 TeV",x_pos_legend,y_pos_legend-0.16*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"|#eta| < 0.8",x_pos_legend,y_pos_legend-0.24*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            if(fCheckBox_sel[2]->GetState() == kButtonDown) plotTopLegend((char*)HistName.Data(),x_pos_legend,y_pos_legend-0.32*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        }
        if(iPad == 1) // phi
        {
            plotTopLegend((char*)"30-40%",x_pos_legend,y_pos_legend-0.08*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"2.76 TeV",x_pos_legend,y_pos_legend-0.16*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"|y| < 0.5",x_pos_legend,y_pos_legend-0.24*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            if(fCheckBox_sel[2]->GetState() == kButtonDown) plotTopLegend((char*)HistName.Data(),x_pos_legend,y_pos_legend-0.32*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

            leg_dNdpT_vs_pT = new TLegend(0.02,0.61,0.58,0.82); // x1,y1,x2,y2
            leg_dNdpT_vs_pT->SetBorderSize(0);
            leg_dNdpT_vs_pT->SetFillColor(0);
            leg_dNdpT_vs_pT->SetTextSize(legend_text_size*scaling_factor_3X3);
            leg_dNdpT_vs_pT->AddEntry((TGraphAsymmErrors*)tg_dNdpT_BW_ana_pid_plot_range[i_mass]->Clone(),"fit","l");
            leg_dNdpT_vs_pT->AddEntry((TGraphAsymmErrors*)tg_dNdpT_BW_ana_pid_plot[i_mass]->Clone(),"prediction","l");
            leg_dNdpT_vs_pT->Draw();
        }
        if(iPad == 2) // d
        {
            plotTopLegend((char*)"20-40%",x_pos_legend,y_pos_legend-0.08*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"2.76 TeV",x_pos_legend,y_pos_legend-0.16*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"|y| < 0.5",x_pos_legend,y_pos_legend-0.24*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            if(fCheckBox_sel[2]->GetState() == kButtonDown) plotTopLegend((char*)HistName.Data(),x_pos_legend,y_pos_legend-0.32*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        }
        if(iPad == 3) // K
        {
            plotTopLegend((char*)"30-40%",x_pos_legend,y_pos_legend-0.08*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"2.76 TeV",x_pos_legend,y_pos_legend-0.16*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"|#eta| < 0.8",x_pos_legend,y_pos_legend-0.24*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            if(fCheckBox_sel[2]->GetState() == kButtonDown) plotTopLegend((char*)HistName.Data(),x_pos_legend,y_pos_legend-0.32*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        }
        if(iPad == 4) // Omega
        {
            plotTopLegend((char*)"20-40%",x_pos_legend,y_pos_legend-0.08*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"2.76 TeV",x_pos_legend,y_pos_legend-0.16*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"|y| < 0.5",x_pos_legend,y_pos_legend-0.24*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            if(fCheckBox_sel[2]->GetState() == kButtonDown) plotTopLegend((char*)HistName.Data(),x_pos_legend,y_pos_legend-0.32*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        }
        if(iPad == 5) // J/Psi
        {
            //plotTopLegend((char*)"20-40%",x_pos_legend,y_pos_legend-0.08*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            //plotTopLegend((char*)"2.76 TeV",x_pos_legend,y_pos_legend-0.16*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            //plotTopLegend((char*)"2.5<y<4",x_pos_legend,y_pos_legend-0.24*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            if(fCheckBox_sel[2]->GetState() == kButtonDown) plotTopLegend((char*)HistName.Data(),x_pos_legend,y_pos_legend-0.32*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        }
        if(iPad == 6) // p
        {
            plotTopLegend((char*)"30-40%",x_pos_legend,y_pos_legend-0.08*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"2.76 TeV",x_pos_legend,y_pos_legend-0.16*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"|#eta| < 0.8",x_pos_legend,y_pos_legend-0.24*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            if(fCheckBox_sel[2]->GetState() == kButtonDown) plotTopLegend((char*)HistName.Data(),x_pos_legend,y_pos_legend-0.32*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        }
        if(iPad == 7) // D0
        {
            plotTopLegend((char*)"30-50%",x_pos_legend,y_pos_legend-0.07*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"2.76 TeV",x_pos_legend,y_pos_legend-0.14*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"|y| < 0.5",x_pos_legend,y_pos_legend-0.21*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            if(fCheckBox_sel[2]->GetState() == kButtonDown) plotTopLegend((char*)HistName.Data(),x_pos_legend,y_pos_legend-0.32*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        }
        if(iPad == 8) // Upsilon
        {
            plotTopLegend((char*)"0-100%",x_pos_legend,y_pos_legend-0.07*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kRed,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"2.76 TeV",x_pos_legend,y_pos_legend-0.14*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"|y|<2.4",x_pos_legend,y_pos_legend-0.21*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            if(fCheckBox_sel[2]->GetState() == kButtonDown) plotTopLegend((char*)HistName.Data(),x_pos_legend,y_pos_legend-0.32*scaling_factor_3X3,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,1,32); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1


            plotTopLegend((char*)"0-100%",2.42,0.375,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,0,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            plotTopLegend((char*)"30-40%",9.24,0.275,legend_text_size*scaling_factor_3X3,kBlack,0.0,42,0,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

            if(!arrow1) arrow1 = new TArrow(4.0,0.35,3.7,0.17,0.01,"|>"); // x1,y1,x2,y2
            arrow1->SetAngle(45.0);
            arrow1->SetLineWidth(3);
            arrow1->SetLineColor(kBlack);
            arrow1->SetFillColor(kBlack);
            arrow1->Draw();

            if(!arrow2) arrow2 = new TArrow(10.0,0.25,9.6,0.18,0.01,"|>"); // x1,y1,x2,y2
            arrow2->SetAngle(45.0);
            arrow2->SetLineWidth(3);
            arrow2->SetLineColor(kBlack);
            arrow2->SetFillColor(kBlack);
            arrow2->Draw();
        }



    }

    c_3X3 ->Modified();
    c_3X3 ->Update();
    printf("dN/dpT plotted \n");

#if 0
    // Check single particle centrality dependence
    if(!c_single_pid)
    {
        c_single_pid = new TCanvas("c_single_pid","c_single_pid",500,20,950,900);
        c_single_pid->SetTopMargin(0.15);
        c_single_pid->SetBottomMargin(0.25);
        c_single_pid->SetRightMargin(0.15);
        c_single_pid->SetLeftMargin(0.25);
        c_single_pid->SetLogy(0);
        c_single_pid->SetTicks(1,1);
        c_single_pid->SetGrid(0,0);
        c_single_pid->SetFillColor(10);
    }

    if(h_frame_single_pid) delete h_frame_single_pid;


    HistName = "h_frame_single_pid";
    h_frame_single_pid = c_single_pid->cd()->DrawFrame(-0.5,-0.1,5.0,2.2,HistName.Data());
    h_frame_single_pid->SetStats(0);
    h_frame_single_pid->SetTitle("");
    h_frame_single_pid->SetTickLength(0.04,"X");
    h_frame_single_pid->SetTickLength(0.04,"Y");
    h_frame_single_pid->GetXaxis()->SetTitleOffset(0.9);
    h_frame_single_pid->GetYaxis()->SetTitleOffset(1.2);
    h_frame_single_pid->GetXaxis()->SetLabelOffset(0.0);
    h_frame_single_pid->GetYaxis()->SetLabelOffset(0.01);
    h_frame_single_pid->GetXaxis()->SetLabelSize(0.08);
    h_frame_single_pid->GetYaxis()->SetLabelSize(0.08);
    h_frame_single_pid->GetXaxis()->SetTitleSize(0.08);
    h_frame_single_pid->GetYaxis()->SetTitleSize(0.08);
    h_frame_single_pid->GetXaxis()->SetNdivisions(505,'N');
    h_frame_single_pid->GetYaxis()->SetNdivisions(505,'N');
    h_frame_single_pid->GetXaxis()->CenterTitle();
    h_frame_single_pid->GetYaxis()->CenterTitle();
    h_frame_single_pid->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_frame_single_pid->GetYaxis()->SetTitle("dN/dp_{T} (GeV/c)^{-1}");
    //h_frame_single_pid->GetYaxis()->SetRangeUser(-0.1,2.2);


    for(Int_t i = 0; i < 6; i++)
    {
        vec_tgae_pT_spectra[1][i] ->SetMarkerColor(i+1);
        vec_tgae_pT_spectra[1][i] ->SetMarkerStyle(20);
        vec_tgae_pT_spectra[1][i] ->SetMarkerSize(0.5);
        vec_tgae_pT_spectra[1][i] ->Draw("same P");
    }

    c_single_pid ->Modified();
    c_single_pid ->Update();
    printf("dN/dpT single plotted \n");
#endif
}
//______________________________________________________________________________



//______________________________________________________________________________
void TBlastWaveGUI::CalcMaxPtLimits()
{
    Double_t beta_max = NEntry_set_limits->GetNumberEntry()->GetNumber();
    if(beta_max < 0.0 || beta_max >= 1.0)
    {
        printf("WARNING: Wrong beta value assigned, set to 0.5");
        beta_max = 0.5;
    }
    calFitRange(beta_max);
    for(int i_pid = 0; i_pid < N_masses_all; ++i_pid)
    {
        if (i_pid < 11 ) arr_NEntry_limits[1][i_pid] ->SetNumber(pT_fit_max[i_pid]);
        if (i_pid >= 11) arr_NEntry_limits_A[1][i_pid] ->SetNumber(pT_fit_max[i_pid]);
    }
    DoSlider();
}
//______________________________________________________________________________



//______________________________________________________________________________
void TBlastWaveGUI::DoSave()
{
    printf("TBlastWaveGUI::DoSave() \n");

    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_save->ChangeBackground(green);

    outputfile->cd();
    gr1 ->Write();
    gr2 ->Write();
    if(c_3X1)
    {
        c_3X1 ->Write();
        c_3X1 ->SaveAs("v2_vs_pT_BW_fits.png");
        c_3X3 ->Write();
        c_3X3 ->SaveAs("dNdpT_vs_pT_BW_fits.png");
    }

}
//______________________________________________________________________________



//______________________________________________________________________________
void TBlastWaveGUI::WriteParams()
{
    printf("WriteParams \n");

    

    //cout  <<"temp " << fit_params[0]  <<" rho0 " << fit_params[1]  << " rho2 " << fit_params[2]  <<" RxOverRy " << fit_params[3] <<endl;
    for (Int_t i_mass =0; i_mass < N_masses_all; i_mass++)
    {
        h_min_val_pT -> AddBinContent(i_mass+1, vec_min_val_pT[i_mass]);
        h_max_val_pT -> AddBinContent(i_mass+1, vec_max_val_pT[i_mass]);
    }
    for (Int_t i_params=0; i_params < 4; i_params++)
    {
        h_fit_params -> AddBinContent(i_params+1, fit_params[i_params]);
    }

    Double_t temp =  h_fit_params -> GetBinContent(1);
    Double_t rho =  h_fit_params -> GetBinContent(2);
    Double_t rho_2 =  h_fit_params -> GetBinContent(3);
    Double_t R =  h_fit_params -> GetBinContent(4);
    cout  <<"T: " <<temp  <<", rho0: " <<rho  << ", rho2: " << rho_2  <<", RxOverRy: " << R  <<endl;

    Double_t min_pt_Pi = h_min_val_pT -> GetBinContent(1);
    Double_t min_pt_t = h_min_val_pT -> GetBinContent(22);
    Double_t max_pt_Pi = h_max_val_pT -> GetBinContent(1);
    Double_t max_pt_t = h_max_val_pT -> GetBinContent(22);
    cout  <<"min_pt_pi+: " << min_pt_Pi  <<", min_pt_t: " << min_pt_t  << ", max_pt_pi+: " << max_pt_Pi  <<", max_pt_t: " << max_pt_t  <<endl;

    TString output_file_name;
    //output_file_name.Clear();
    output_file_name = combo_energy.Data();
    output_file_name += "_";
    output_file_name += cent_lower.Data();
    output_file_name += "_";
    output_file_name += cent_upper.Data();
    for ( Int_t i_mass = 0; i_mass < N_masses_all; i_mass ++)
    {
        if (!fCheckBox_pid[i_mass] ->IsDown()) continue;
        output_file_name += "_";
        if (i_mass == 15)
        {
            if (fCheckBox_pid[i_mass] ->IsDown() || fCheckBox_pid_fit_dNdpt[i_mass] ->IsDown()) output_file_name += "JPsi";

        }
        if (i_mass != 15)
        {
            if (fCheckBox_pid[i_mass] ->IsDown() || fCheckBox_pid_fit_dNdpt[i_mass] ->IsDown()) output_file_name += label_full_pid_spectra[i_mass];
        }
    }
    output_file_name += "_";
    output_file_name += Add_output_file_name->GetDisplayText();
    output_file_name += ".root";
    cout << "output_file_name: " << output_file_name.Data() << endl;
    RootFileFitParams = new TFile(output_file_name.Data(), "RECREATE");
    if ( RootFileFitParams->IsOpen() ) printf("File opened successfully\n");
    RootFileFitParams->cd();
    h_fit_params->Write("h_fit_params");
    h_min_val_pT->Write("h_min_val_pT");
    h_max_val_pT->Write("h_max_val_pT");


}
//______________________________________________________________________________


void v2_slider()
{
    new TBlastWaveGUI();
}
