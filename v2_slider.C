
#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"
#include "TGLayout.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGTextEntry.h"
#include "TGTripleSlider.h"
#include <TGSlider.h>

#include "./functions_BW.h"

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

    TProfile* tp_v2_vs_pT_mesons[8][9][9][9][9][9]; // [i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]
    TH1F*     h_dN_dpT_mesons[8][9][9][9][9][9]; // [i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]
    TList* arr_list[9][9];
    vector<TGHorizontalFrame*>  vec_Hframe;
    vector<TGHSlider*>          vec_slider;
    vector<TGLayoutHints*>      vec_LayoutHints;
    vector<TGTextEntry*>        vec_TextEntry;
    vector<TGTextBuffer*>       vec_TextBuffer;
    vector<TGLabel*>            vec_TGLabel;

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
    TGVerticalFrame *hVframeD5a[4];
    TGVerticalFrame *hVframeD3;
    TGVerticalFrame *hVframeD4;
    TGVerticalFrame *vframeD1;
    TGCompositeFrame *cframe2;
    TGTextButton *ButtonD1a;
    TGTextButton *Button_exit;
    TGTextButton *Button_save;
    TGHProgressBar* fHProg1 = NULL;
    TGHProgressBar* fHProg2 = NULL;
    TGVProgressBar* fVProg1;
    TGLayoutHints* fHint2;
    TGTextButton      *fGO;
    TGLayoutHints* fHint3;
    TGLayoutHints* LHintsD4a;
    TGNumberEntry* NEntryD3a;
    TGLabel*       LabelD4a;
    TGNumberEntry* NEntryD3b;
    TGLabel*       LabelD4b;
    TGTransientFrame* frame_TGTransient;
    TGTransientFrame* frame_TGTransientB;
    TGTextButton      *Button_take_params_MC_to_ana;
    TGTextButton      *Button_take_params_Set_to_ana;

    TGCheckButton* fCheckBox_sel[2];
    TGCheckButton* fCheckBox_pid[8];
    TGCheckButton* fCheckBox_v2_dNdpT[2];
    TGLayoutHints* fLCheckBox;
    TString label_checkbox[2] = {"Plot MC","Plot Ana"};

    TGHorizontalFrame* arr_HFrame_NEntry_limits[8];
    TGVerticalFrame* arr_VFrame_NEntry_limits[2];
    TGLabel*           arr_Label_NEntry_limits[2][8];
    TGNumberEntry*     arr_NEntry_limits[2][8];
    TGNumberEntry*     arr_NEntry_ana_params[4];
    TGLabel*           arr_Label_NEntry_ana_params[4];
    TGLayoutHints* arr_fL1[2];
    Double_t min_max_pT_range_pid[2][8];
    TGraph* tg_v2_BW_ana_pid[8];
    TGraph* tg_dNdpT_BW_ana_pid[8];

    TGraph* tg_v2_BW_ana_pid_min[8];
    TGraph* tg_dNdpT_BW_ana_pid_min[8];

    TGraph* tg_v2_BW_ana_pid_plot[8];
    TGraph* tg_v2_BW_ana_pid_plot_range[8];
    TGraph* tg_dNdpT_BW_ana_pid_plot[8];


    Double_t var_test = 5.2;

    TGTextButton *Button_minimize;
    TGTextButton *Button_stop_minimize;

    TGTextButton *Button_minimize_ana;
    TGTextButton *Button_stop_minimize_ana;

    TGTextButton *Button_make_plot_v2;
    TGTextButton *Button_make_plot_dNdpT;

    TFile* inputfile;
    TH1D* h_dummy;
    TH1D* h_dummy_plot_v2[2] = {NULL,NULL};
    TH1D* h_dummy_dNdpT;
    TLegend* leg_v2_vs_pT_A = NULL;
    TLegend* leg_v2_vs_pT_B = NULL;
    Int_t flag_stop_minimize = 0;
    Int_t flag_minimization_ana = 0;
    Double_t integration_range_pid[8][2] = {0.0};

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
    TLatex* TLatex_legend_dNdpT[8];
    TLatex* TLatex_legend_v2_plot[8];

    TLine* TL_line_base = NULL;
    TLine* TL_line_base_plot[2] = {NULL,NULL};

    Double_t T_BW_fit_ana        = -1.0;
    Double_t Rho0_BW_fit_ana     = -1.0;
    Double_t Rho2_BW_fit_ana     = -1.0;
    Double_t RxOverRy_BW_fit_ana = -1.0;

    TH1F* h_frame_2X1[2] = {NULL,NULL};
    TCanvas* c_2X1 = NULL;
    TGraph* tg_label_plot[3][8];

    TFile* outputfile;

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
    void DoSave();
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
    min_max_pT_range_pid[0][0] = 0.1; // pi
    min_max_pT_range_pid[1][0] = 1.6;
    min_max_pT_range_pid[0][1] = 0.1; // K
    min_max_pT_range_pid[1][1] = 1.9;
    min_max_pT_range_pid[0][2] = 0.1; // p
    min_max_pT_range_pid[1][2] = 2.2;
    min_max_pT_range_pid[0][3] = 0.1; // phi
    min_max_pT_range_pid[1][3] = 2.2;
    min_max_pT_range_pid[0][4] = 0.1; // Omega
    min_max_pT_range_pid[1][4] = 4.0;
    min_max_pT_range_pid[0][5] = 0.1; // D0
    min_max_pT_range_pid[1][5] = 8.5;
    min_max_pT_range_pid[0][6] = 0.1; // J/Psi
    min_max_pT_range_pid[1][6] = 15.0;
    min_max_pT_range_pid[0][7] = 0.1; // Upsilon
    min_max_pT_range_pid[1][7] = 15.0;
    //--------------------------------------------------------------------


    for(int i_mass = 0; i_mass < N_masses; ++i_mass)
    {
        tg_v2_BW_ana_pid[i_mass]    = NULL;
        tg_dNdpT_BW_ana_pid[i_mass] = NULL;

        tg_v2_BW_ana_pid_min[i_mass]    = NULL;
        tg_dNdpT_BW_ana_pid_min[i_mass] = NULL;

        tg_v2_BW_ana_pid_plot[i_mass]    = NULL;
        tg_v2_BW_ana_pid_plot_range[i_mass]    = NULL;
        tg_dNdpT_BW_ana_pid_plot[i_mass] = NULL;

        TLatex_legend_dNdpT[i_mass] = NULL;
        TLatex_legend_v2_plot[i_mass] = NULL;

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
    make_5_60_spectra();
    Init_v2_Mathematica();


    inputfile = TFile::Open("./Data/merge_out_v2_boost.root");

    for(Int_t i_R_x = 0; i_R_x < 9; i_R_x++)
    {
        for(Int_t i_fboost = 0; i_fboost < 9; i_fboost++)
        {
            printf("i_R_x: %d, i_fboost: %d \n",i_R_x,i_fboost);
            arr_list[i_R_x][i_fboost] = (TList*)inputfile->Get(Form("list_BW_Rx%d_fb%d",i_R_x,i_fboost));
        }
    }


    // Determine integration range of dNdpT for data
    for(Int_t i_mass = 0; i_mass < 8; i_mass++)
    {
        Double_t x_val_data_first, y_val_data_first, x_val_data_last, y_val_data_last, x_err_low_data, x_err_high_data;
        tgae_dN_dpT_mesons_data[i_mass] ->GetPoint(0,x_val_data_first,y_val_data_first);
        x_err_low_data = tgae_dN_dpT_mesons_data[i_mass] ->GetErrorXlow(0);
        tgae_dN_dpT_mesons_data[i_mass] ->GetPoint(tgae_dN_dpT_mesons_data[i_mass] ->GetN()-1,x_val_data_last,y_val_data_last);
        x_err_high_data = tgae_dN_dpT_mesons_data[i_mass] ->GetErrorXhigh(tgae_dN_dpT_mesons_data[i_mass] ->GetN()-1);

        integration_range_pid[i_mass][0] = x_val_data_first - x_err_low_data;
        integration_range_pid[i_mass][1] = x_val_data_last  + x_err_high_data;
    }


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
                        for(Int_t i_mass = 0; i_mass < 8; i_mass++)
                        {
                            Int_t i_N = i_mass + 8*i_rho_a + 8*9*i_rho_0 + 8*9*9*i_Temp;
                            Int_t i_pos_v2    = i_N*2;
                            Int_t i_pos_dNdpT = i_N*2 + 1;
                            tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] = (TProfile*)arr_list[i_R_x][i_fboost] ->At(i_pos_v2);
                            h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a]    = (TH1F*)arr_list[i_R_x][i_fboost]     ->At(i_pos_dNdpT);

                            // integrate only in range of available data
                            Int_t start_integrate_bin = h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->FindBin(integration_range_pid[i_mass][0]);
                            Int_t stop_integrate_bin  = h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->FindBin(integration_range_pid[i_mass][1]);
                            Double_t integral = h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->Integral(start_integrate_bin,stop_integrate_bin,"width");
                            if(integral <= 0.0) continue;
                            h_dN_dpT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->Scale(1.0/integral);
                        }
                    }
                }
            }
        }
    }

    char buf[32];
    SetCleanup(kDeepCleanup);
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

    //----------------------------------------------------
    // Slider start XA
    //Int_t start_pos_slider[5] = {5,7,3,5,0};
    //Int_t start_pos_slider[5] = {5,4,3,5,0};
    //Int_t start_pos_slider[5] = {4,4,2,6,0}; // overall OK
    //Int_t start_pos_slider[5] = {4,7,3,6,0}; // very good for v2 only
    Int_t start_pos_slider[5] = {3,5,3,7,0}; // overall quite good

    // Add a slider
    for(Int_t i_param = 0; i_param < 5; i_param++)
    {
        vec_Hframe[i_param] = new TGHorizontalFrame(this, 0, 0);
        vec_slider[i_param] = new TGHSlider(vec_Hframe[i_param],250,kSlider1|kScaleDownRight,1);
        vec_slider[i_param]->Connect("PositionChanged(Int_t)", "TBlastWaveGUI",this, "DoSlider()");
        vec_slider[i_param]->SetRange(0,8);
        vec_slider[i_param]->SetPosition(start_pos_slider[i_param]);
        vec_LayoutHints[i_param] = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 1, 1, 1, 1); // handles size of slider and number box
        vec_Hframe[i_param]->AddFrame(vec_slider[i_param], vec_LayoutHints[i_param]);
        AddFrame(vec_Hframe[i_param], vec_LayoutHints[i_param]);


        // Add number text box
        vec_TextEntry[i_param] = new TGTextEntry(vec_Hframe[i_param], vec_TextBuffer[i_param] = new TGTextBuffer(5), HId1);
        //vec_TextEntry[i_param] ->SetDefaultSize(0,0);
        vec_TextEntry[i_param] ->SetTextColor(kRed);
        vec_TextEntry[i_param]->SetToolTipText("Minimum (left) Value of Slider");
        vec_TextBuffer[i_param]->AddText(0, "0.0");
        vec_TextEntry[i_param]->Connect("TextChanged(char*)", "TBlastWaveGUI", this,"DoText(char*)");
        //vec_Hframe[i_param]->Resize(100, 25);
        vec_Hframe[i_param]->AddFrame(vec_TextEntry[i_param], vec_LayoutHints[i_param]);
        AddFrame(vec_Hframe[i_param], vec_LayoutHints[i_param]);
        vec_TGLabel[i_param] = new TGLabel(vec_Hframe[i_param], vec_label[i_param].Data());
        vec_Hframe[i_param]->AddFrame(vec_TGLabel[i_param], new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

        vec_Hframe[i_param] ->SetHeight(1);
        //vec_Hframe[i_param] ->DrawBorder();
    }
    //----------------------------------------------------


    SetWindowName("Elliptic flow data");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();
    this ->Resize(1100,1100);
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
    fCanvasB->GetCanvas()->Divide(4,2);
    for(Int_t iPad = 1; iPad <= 8; iPad++)
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

    for(Int_t iPad = 1; iPad <= 8; iPad++)
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
    hframeD1  = new TGHorizontalFrame(FrameD,200,100);

    // exit button
    Button_exit = new TGTextButton(hframeD1, "&Exit ","gApplication->Terminate(0)");
    hframeD1->AddFrame(Button_exit, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // save button
    Button_save = new TGTextButton(hframeD1, "&Save ",10);
    Button_save->Connect("Clicked()", "TBlastWaveGUI", this, "DoSave()");
    hframeD1->AddFrame(Button_save, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    FrameD ->AddFrame(hframeD1, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



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



    //--------------
    // A horizontal frame
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

    // exit button
    Button_minimize_ana = new TGTextButton(hframeD2b, "Minimize Ana",10);
    Button_minimize_ana->Connect("Clicked()", "TBlastWaveGUI", this, "DoMinimize_ana()");

    Button_minimize_ana->ChangeBackground(red);

    hframeD2b->AddFrame(Button_minimize_ana, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // save button
    Button_stop_minimize_ana = new TGTextButton(hframeD2b, "Stop minimize Ana",10);
    Button_stop_minimize_ana->Connect("Clicked()", "TBlastWaveGUI", this, "StopMinimize()");
    hframeD2b->AddFrame(Button_stop_minimize_ana, new TGLayoutHints(kLHintsCenterX,5,5,3,4));


    FrameD ->AddFrame(hframeD2b, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    //--------------
    hframeD3  = new TGHorizontalFrame(FrameD,200,100);
    fLCheckBox = new TGLayoutHints(kLHintsTop | kLHintsLeft,0, 0, 5, 0);
    for(Int_t i_particle = 0; i_particle < 8; i_particle++)
    {
        fCheckBox_pid[i_particle]  = new TGCheckButton(hframeD3, new TGHotString(label_full_pid_spectra[i_particle].Data()), -1);
        fCheckBox_pid[i_particle] ->SetState(kButtonDown);
        fCheckBox_pid[i_particle] ->Connect("Clicked()", "TBlastWaveGUI", this, "DoSlider()");
        hframeD3->AddFrame(fCheckBox_pid[i_particle], fLCheckBox);
    }
    FrameD ->AddFrame(hframeD3, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    //--------------
    hframeD4  = new TGHorizontalFrame(FrameD,200,100);
    for(Int_t i_v2_dNdpT = 0; i_v2_dNdpT < 2; i_v2_dNdpT++)
    {
        fCheckBox_v2_dNdpT[i_v2_dNdpT]  = new TGCheckButton(hframeD4, new TGHotString(label_v2_dNdpT[i_v2_dNdpT].Data()), -1);
        fCheckBox_v2_dNdpT[i_v2_dNdpT] ->SetState(kButtonDown);
        hframeD4->AddFrame(fCheckBox_v2_dNdpT[i_v2_dNdpT], fLCheckBox);
    }
    for(Int_t i_cb = 0; i_cb < 2; i_cb++)
    {
        fCheckBox_sel[i_cb]  = new TGCheckButton(hframeD4, new TGHotString(label_checkbox[i_cb].Data()), -1);
        fCheckBox_sel[i_cb] ->SetState(kButtonDown);
        fCheckBox_sel[i_cb] ->Connect("Clicked()", "TBlastWaveGUI", this, "DoSlider()");
        hframeD4->AddFrame(fCheckBox_sel[i_cb], fLCheckBox);
    }
    FrameD ->AddFrame(hframeD4, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    //--------------
    TString arr_label_params_ana[4] = {"SET T","SET rho0","SET rhoa","SET Rx"};
    hframeD5a  = new TGHorizontalFrame(FrameD,200,100);
    for(Int_t i_param = 0; i_param < 4; i_param++)
    {
        hVframeD5a[i_param] = new TGVerticalFrame(hframeD5a, 200,200);
        arr_NEntry_ana_params[i_param] = new TGNumberEntry(hVframeD5a[i_param], 0.0, 12,(TGNumberFormat::EStyle) 2);
        arr_NEntry_ana_params[i_param] ->SetNumStyle( TGNumberFormat::kNESRealTwo); // https://root.cern.ch/doc/master/classTGNumberFormat.html#a8a0f81aac8ac12d0461aef554c6271ad
        hVframeD5a[i_param]->AddFrame(arr_NEntry_ana_params[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));

        TString label_entry = arr_label_params_ana[i_param];
        arr_Label_NEntry_ana_params[i_param] = new TGLabel(hVframeD5a[i_param], label_entry.Data(), myGC(), myfont->GetFontStruct());
        hVframeD5a[i_param]->AddFrame(arr_Label_NEntry_ana_params[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
        hframeD5a->AddFrame(hVframeD5a[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    }
    FrameD ->AddFrame(hframeD5a, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    //--------------
    hframeD5  = new TGHorizontalFrame(FrameD,200,100);

    Button_take_params_MC_to_ana = new TGTextButton(hframeD5, "Use MC params for Ana",10);
    Button_take_params_MC_to_ana->Connect("Clicked()", "TBlastWaveGUI", this, "TakeParamsFromMC()");
    hframeD5->AddFrame(Button_take_params_MC_to_ana, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    Button_take_params_Set_to_ana = new TGTextButton(hframeD5, "Use SET params for Ana",10);
    Button_take_params_Set_to_ana->Connect("Clicked()", "TBlastWaveGUI", this, "TakeParamsFromSet()");
    hframeD5->AddFrame(Button_take_params_Set_to_ana, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    FrameD ->AddFrame(hframeD5, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------




    //--------------
    hframeD6  = new TGHorizontalFrame(FrameD,200,100);

    Button_make_plot_v2 = new TGTextButton(hframeD6, "Make plot v2",10);
    Button_make_plot_v2->Connect("Clicked()", "TBlastWaveGUI", this, "MakePlotv2()");
    hframeD6->AddFrame(Button_make_plot_v2, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    Button_make_plot_dNdpT = new TGTextButton(hframeD6, "Make plot dN/dpT",10);
    Button_make_plot_dNdpT->Connect("Clicked()", "TBlastWaveGUI", this, "MakePlotdNdpT()");
    hframeD6->AddFrame(Button_make_plot_dNdpT, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    FrameD ->AddFrame(hframeD6, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------





    FrameD ->MapSubwindows();
    FrameD ->MapWindow();
    FrameD ->Resize(450,400); // size of frame
    FrameD ->Move(1250,750); // position of frame



    //--------------
    LHintsD4a = new TGLayoutHints(kLHintsCenterX,5,5,3,4);

    TString arr_label_pid[8]     = {"pi","K","p","phi","Omega","D0","J/Psi","Y"};
    TString arr_label_min_max[2] = {"min","max"};

    arr_fL1[0] = new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2);
    arr_fL1[1] = new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2);

    Double_t max_val_VF[2] = {300,500};

    frame_TGTransient = new TGTransientFrame(gClient->GetRoot(), FrameD, 10, 10, kHorizontalFrame);
    FrameD ->AddFrame(frame_TGTransient,LHintsD4a);
    for(Int_t i_min_max = 0; i_min_max < 2; i_min_max++)
    {
        //arr_VFrame_NEntry_limits[i_min_max] = new TGVerticalFrame(FrameD, 200, max_val_VF[i_min_max]);
        //FrameD->AddFrame(arr_VFrame_NEntry_limits[i_min_max], arr_fL1[i_min_max]);

        arr_VFrame_NEntry_limits[i_min_max] = new TGVerticalFrame(frame_TGTransient, 200, max_val_VF[i_min_max]);
        frame_TGTransient->AddFrame(arr_VFrame_NEntry_limits[i_min_max], arr_fL1[i_min_max]);

        for(Int_t i_pid = 0; i_pid < 8; i_pid++)
        {
            arr_NEntry_limits[i_min_max][i_pid] = new TGNumberEntry(arr_VFrame_NEntry_limits[i_min_max], min_max_pT_range_pid[i_min_max][i_pid], 12,(TGNumberFormat::EStyle) 1);
            arr_NEntry_limits[i_min_max][i_pid] ->Connect("ValueSet(Long_t)", "TBlastWaveGUI", this, "DoSlider()");
            arr_NEntry_limits[i_min_max][i_pid]->SetNumStyle( TGNumberFormat::kNESRealOne); // https://root.cern.ch/doc/master/classTGNumberFormat.html#a8a0f81aac8ac12d0461aef554c6271ad
            arr_VFrame_NEntry_limits[i_min_max]->AddFrame(arr_NEntry_limits[i_min_max][i_pid], LHintsD4a);
            TString label_entry = "pT " + arr_label_min_max[i_min_max] + " " + arr_label_pid[i_pid];
            arr_Label_NEntry_limits[i_min_max][i_pid] = new TGLabel(arr_VFrame_NEntry_limits[i_min_max], label_entry.Data(), myGC(), myfont->GetFontStruct());
            arr_VFrame_NEntry_limits[i_min_max]->AddFrame(arr_Label_NEntry_limits[i_min_max][i_pid], LHintsD4a);
        }
    }

    frame_TGTransient->MapSubwindows();
    frame_TGTransient->Resize();
    frame_TGTransient->CenterOnParent();
    frame_TGTransient->SetWindowName("pT limits");
    frame_TGTransient->MapWindow();
    frame_TGTransient ->Move(1250,250); // position of frame
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
    //cout << "DoSlider started" << endl;

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

    for(Int_t i_mass = 0; i_mass < 8; i_mass++)
    {
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

    if(TL_line_base) delete TL_line_base;
    TL_line_base = PlotLine(0.0,15.0,0.0,0.0,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)

    for(Int_t i_mass = 0; i_mass < 8; i_mass++)
    {
        if(fCheckBox_pid[i_mass]->GetState() == kButtonDown)
        {
            //printf("Draw i_mass: %d \n", i_mass);
            tgae_v2_vs_pT_mesons_data[i_mass] ->SetMarkerColor(arr_color_mass[i_mass]);
            tgae_v2_vs_pT_mesons_data[i_mass] ->Draw("same P");

            if(fCheckBox_sel[0]->GetState() == kButtonDown)
            {
                tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineWidth(4);
                tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineColor(arr_color_mass[i_mass]);
                tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->Draw("same L hist");
            }
        }
    }

    if(fCheckBox_sel[1]->GetState() == kButtonDown && flag_minimization_ana)
    {
        for(int i_mass = 0; i_mass < N_masses; ++i_mass)
        {
            if(fCheckBox_pid[i_mass]->GetState() == kButtonDown)
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
    leg_v2_vs_pT_B->AddEntry((TGraphAsymmErrors*)tgae_v2_vs_pT_mesons_data[6]->Clone(),"J/#Psi","p");
    leg_v2_vs_pT_B->AddEntry((TGraphAsymmErrors*)tgae_v2_vs_pT_mesons_data[7]->Clone(),"#Upsilon","p");
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
    Double_t x_range_dNdpT[8] = {3.5,4.1,4.5,5.0,5.0,5.0,8.5,15.0};
    Double_t y_range_dNdpT[8] = {2.2,1.4,1.1,1.0,1.0,1.0,0.45,0.35};

    for(Int_t i_mass = 0; i_mass < 8; i_mass++)
    {
        fCanvasB ->GetCanvas()->cd(i_mass + 1);
        h_dummy_dNdpT ->GetXaxis()->SetRangeUser(0.0,x_range_dNdpT[i_mass]);
        h_dummy_dNdpT ->GetYaxis()->SetRangeUser(-0.2,y_range_dNdpT[i_mass]);
        h_dummy_dNdpT ->DrawCopy("");

        tgae_dN_dpT_mesons_data[i_mass] ->SetMarkerColor(kBlack);
        tgae_dN_dpT_mesons_data[i_mass] ->SetMarkerStyle(20);
        tgae_dN_dpT_mesons_data[i_mass] ->SetMarkerSize(1.0);
        tgae_dN_dpT_mesons_data[i_mass] ->Draw("same P");

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
void TBlastWaveGUI::DoMinimize_ana()
{
    cout << "DoMinimize_ana started" << endl;
    N_calls_BW_ana = 0;

    fHProg2 ->Reset();

    Pixel_t yellow;
    gClient->GetColorByName("yellow", yellow);
    Button_minimize_ana->ChangeBackground(yellow);

    Button_take_params_MC_to_ana->ChangeBackground(yellow);


    auto chi2Function = [&](const Double_t *par)
    {
        //minimisation function computing the sum of squares of residuals
        N_calls_BW_ana++;

        double chi2 = 0;

        const double T = par[0];       // fit parameter: Temp in GeV
        const double rho0 = par[1];   // fit parameter: transverse rapidity
        const double rho2 = par[2];    // fit parameter: azimuthal modulation of transverse rapidity
        const double RxOverRy = par[3]; // fit parameter: ratio of the radii Rx and Ry of the freeze-out ellipse in the transverse plane

        for(int i_mass = 0; i_mass < N_masses; ++i_mass)
        {
            if(!(fCheckBox_pid[i_mass]->GetState() == kButtonDown)) continue;
            Double_t min_val_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();
            Double_t max_val_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();

            const double m = arr_quark_mass_meson[i_mass];       // in GeV

            if(tg_v2_BW_ana_pid_min[i_mass])    delete tg_v2_BW_ana_pid_min[i_mass];
            if(tg_dNdpT_BW_ana_pid_min[i_mass]) delete tg_dNdpT_BW_ana_pid_min[i_mass];

            tg_v2_BW_ana_pid_min[i_mass]    = new TGraph();
            tg_dNdpT_BW_ana_pid_min[i_mass] = new TGraph();

            //--------------------------------------------------
            // v2 chi2
            if((fCheckBox_v2_dNdpT[0]->GetState() == kButtonDown))
            {
                Int_t i_point_ana = 0;
                for(int i_point = 0; i_point < tgae_v2_vs_pT_mesons_data[i_mass]->GetN(); ++i_point)
                {
                    double v2_data = 0.0;
                    double pt_data = 0.0;

                    tgae_v2_vs_pT_mesons_data[i_mass]->GetPoint(i_point,pt_data,v2_data);
                    double v2_err = tgae_v2_vs_pT_mesons_data[i_mass]->GetErrorYhigh(i_point);
                    // cout << "i_point = " << i_point << ", pt_data = " << pt_data << "v2 = " << v2_data << " +/- " << v2_err << endl;

                    if(pt_data > max_val_pT) break;

                    if(pt_data > min_val_pT)
                    {
                        double v2_BW = 0;
                        double inv_yield_BW = 0;

                        // blast wave parameters
                        const double pt_BW = pt_data;         // in GeV

                        blastwave_yield_and_v2(pt_BW, m, T, rho0, rho2, RxOverRy, inv_yield_BW, v2_BW);

                        tg_v2_BW_ana_pid_min[i_mass] ->SetPoint(i_point_ana,pt_BW,v2_BW);

                        // cout << "i_point = " << i_point << ", pt_BW = " << pt_BW << ", v2_BW = " << v2_BW << endl;
                        double diff   = (v2_data - v2_BW)/v2_err;
                        //double diff_yield = (inv_yield_BW - inv_yield_data)/yield_err;
                        chi2 += diff*diff;
                        i_point_ana++;
                    }
                }
            }
            //--------------------------------------------------


            //--------------------------------------------------
            // dNdpT chi2
            if((fCheckBox_v2_dNdpT[1]->GetState() == kButtonDown))
            {
                Double_t integral_ana = 0.0;
                vector< vector<Double_t> > vec_data_BW;
                vec_data_BW.resize(4); // data, BW, err, pT
                if(i_mass != 7) // Upsilon spectra do not exist
                {
                    // First the integral of BW needs to be calculated to do a shape comparison to the data
                    for(int i_point = 0; i_point < tgae_dN_dpT_mesons_data[i_mass]->GetN(); ++i_point)
                    {
                        double dNdpT_data = 0.0;
                        double pt_data    = 0.0;

                        tgae_dN_dpT_mesons_data[i_mass]->GetPoint(i_point,pt_data,dNdpT_data);
                        double dNdpT_err = tgae_dN_dpT_mesons_data[i_mass]->GetErrorYhigh(i_point);

                        Double_t x_err_low_data  = tgae_dN_dpT_mesons_data[i_mass] ->GetErrorXlow(i_point);
                        Double_t x_err_high_data = tgae_dN_dpT_mesons_data[i_mass] ->GetErrorXhigh(i_point);
                        Double_t bin_width = x_err_low_data + x_err_high_data;

                        // cout << "i_point = " << i_point << ", pt_data = " << pt_data << "v2 = " << dNdpT_data << " +/- " << dNdpT_err << endl;

                        //if(pt_data > max_val_pT) break;

                        //if(pt_data > min_val_pT)
                        {
                            double v2_BW = 0;
                            double inv_yield_BW = 0;

                            // blast wave parameters
                            const double pt_BW = pt_data;         // in GeV

                            blastwave_yield_and_v2(pt_BW, m, T, rho0, rho2, RxOverRy, inv_yield_BW, v2_BW); // invariant yield: (1/pT) (dN/dpT)
                            vec_data_BW[0].push_back(dNdpT_data);
                            vec_data_BW[1].push_back(inv_yield_BW*pt_BW);
                            vec_data_BW[2].push_back(dNdpT_err);
                            vec_data_BW[3].push_back(pt_BW);

                            integral_ana += inv_yield_BW*bin_width*pt_BW;
                            // cout << "i_point = " << i_point << ", pt_BW = " << pt_BW << ", v2_BW = " << v2_BW << endl;
                            //double diff   = (dNdpT_data - inv_yield_BW)/dNdpT_err;
                            //double diff_yield = (inv_yield_BW - inv_yield_data)/yield_err;
                            //chi2 += diff*diff;
                        }
                    }


                    // Calculate chi2 for dNdpT
                    if(integral_ana > 0.0)
                    {
                        Int_t i_point_ana = 0;
                        for(Int_t i_point = 0; i_point < (Int_t)vec_data_BW[0].size(); i_point++)
                        {
                            Double_t pt_data = vec_data_BW[3][i_point];
                            if(pt_data > max_val_pT) break;

                            if(pt_data > min_val_pT)
                            {
                                // Normalize BW to integral within range of data
                                double diff   = (vec_data_BW[0][i_point] - (vec_data_BW[1][i_point]/integral_ana))/vec_data_BW[2][i_point];
                                chi2 += diff*diff;
                                tg_dNdpT_BW_ana_pid_min[i_mass] ->SetPoint(i_point_ana,pt_data,vec_data_BW[1][i_point]/integral_ana);
                                i_point_ana++;
                            }
                        }
                    }
                }
            }

            //integration_range_pid[i_mass][0] = x_val_data_first - x_err_low_data;
            //integration_range_pid[i_mass][1] = x_val_data_last  + x_err_high_data;
            //--------------------------------------------------
        }

        if(N_calls_BW_ana % 5 == 0) fHProg2->Increment(5);
        //printf("fraction total: %4.2f%%, fraction added: %4.2f%% \n",fraction_total*100.0,fraction_use*100.0);

        fCanvas->GetCanvas()->cd();
        HistName = "#chi_{Ana}^{2}/ndf = ";
        sprintf(NoP,"%4.2f",chi2);
        HistName += NoP;
        //printf("chi2 (ana): %4.3f \n",chi2);
        if(TLatex_chi2_ana) delete TLatex_chi2_ana;
        TLatex_chi2_ana = plotTopLegend((char*)HistName.Data(),0.18,0.80,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        for(int i_mass = 0; i_mass < N_masses; ++i_mass)
        {
            if(!tg_v2_BW_ana_pid_min[i_mass]) continue;
            tg_v2_BW_ana_pid_min[i_mass] ->SetLineColor(arr_color_mass[i_mass]);
            tg_v2_BW_ana_pid_min[i_mass] ->SetLineWidth(3);
            tg_v2_BW_ana_pid_min[i_mass] ->SetLineStyle(9);
            tg_v2_BW_ana_pid_min[i_mass] ->Draw("same L");
        }

        fCanvas->GetCanvas()->Modified();
        fCanvas->GetCanvas()->Update();

        for(int i_mass = 0; i_mass < N_masses; ++i_mass)
        {
            if(!tg_dNdpT_BW_ana_pid_min[i_mass]) continue;
            fCanvasB ->GetCanvas()->cd(i_mass+1);
            tg_dNdpT_BW_ana_pid_min[i_mass] ->SetLineColor(kAzure-2);
            tg_dNdpT_BW_ana_pid_min[i_mass] ->SetLineWidth(3);
            tg_dNdpT_BW_ana_pid_min[i_mass] ->SetLineStyle(1);
            tg_dNdpT_BW_ana_pid_min[i_mass] ->Draw("same L");
        }

        fCanvasB->GetCanvas()->Modified();
        fCanvasB->GetCanvas()->Update();

        gSystem->Sleep(100);
        gSystem->ProcessEvents();

        if(chi2 < chi2_final_min_ana) chi2_final_min_ana = chi2;

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

    const ROOT::Fit::FitResult & result = fitter.Result();
    result.Print(std::cout);
    const double *fitpar = result.GetParams();
    double T_BW        = fitpar[0];
    double Rho0_BW     = fitpar[1];
    double Rho2_BW     = fitpar[2];
    double RxOverRy_BW = fitpar[3];

    T_BW_fit_ana        = T_BW;
    Rho0_BW_fit_ana     = Rho0_BW;
    Rho2_BW_fit_ana     = Rho2_BW;
    RxOverRy_BW_fit_ana = RxOverRy_BW;


    double Chi2_ana = result.Chi2();
    double NDF_ana  = result.Ndf();
    cout << "T_BW = " << T_BW << ", Rho0_BW = " << Rho0_BW << ", Rho2_BW = " << Rho2_BW << ", RxOverRy_BW = " << RxOverRy_BW << endl;

    Plot_curves_ana(T_BW,Rho0_BW,Rho2_BW,RxOverRy_BW);

    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_minimize_ana->ChangeBackground(green);
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
        for(Int_t i_pid = 0; i_pid < 8; i_pid++)
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


                        for(Int_t i_mass = 0; i_mass < 8; i_mass++) //

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

                            if(i_mass != 7) // Upsilon spectra do not exist
                            {
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


    //--------------------------------------------------------------------
    Double_t min_max_pT_range_pid_plot_ana[2][8];
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
    //--------------------------------------------------------------------


    //------------------------------------------------------
    // v2 plots
    fCanvas ->GetCanvas() ->cd();

    const Int_t N_points_BW_ana = 35;
    for(int i_mass = 0; i_mass < N_masses; ++i_mass)
    {
        Double_t delta_pT = (Double_t)((min_max_pT_range_pid_plot_ana[1][i_mass] - min_max_pT_range_pid_plot_ana[0][i_mass])/(Double_t)N_points_BW_ana);
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
                blastwave_yield_and_v2(pt_BW, arr_quark_mass_meson[i_mass], T_BW, Rho0_BW, Rho2_BW, RxOverRy_BW, inv_yield_BW, v2_BW);
                tg_v2_BW_ana_pid[i_mass] ->SetPoint(i_pT,pt_BW,v2_BW);
            }
        }
        tg_v2_BW_ana_pid[i_mass] -> SetLineColor(arr_color_mass[i_mass]);
        tg_v2_BW_ana_pid[i_mass] -> SetLineWidth(3);
        tg_v2_BW_ana_pid[i_mass] -> SetLineStyle(9);
        //if(!(fCheckBox_pid[i_mass]->GetState() == kButtonDown)) tg_v2_BW_ana_pid[i_mass] -> SetLineStyle(9);
        if(fCheckBox_sel[1]->GetState() == kButtonDown) tg_v2_BW_ana_pid[i_mass] -> Draw("same L");
    }

    printf("v2 ana plotted \n");
    fCanvas->GetCanvas()->Modified();
    fCanvas->GetCanvas()->Update();
    //------------------------------------------------------



    //------------------------------------------------------
    // dNdpT plots
    fCanvasB ->GetCanvas() ->cd();

    for(int i_mass = 0; i_mass < N_masses; ++i_mass)
    {
        //Double_t delta_pT = (Double_t)((integration_range_pid[i_mass][1] - 0.0)/(Double_t)N_points_BW_ana);
        Double_t delta_pT = (Double_t)((min_max_pT_range_pid_plot_ana[1][i_mass] - min_max_pT_range_pid_plot_ana[0][i_mass])/(Double_t)N_points_BW_ana);
        //Double_t min_pT = arr_NEntry_limits[0][i_mass]->GetNumberEntry()->GetNumber();
        //Double_t max_pT = arr_NEntry_limits[1][i_mass]->GetNumberEntry()->GetNumber();

        fCanvasB ->GetCanvas()->cd(i_mass + 1);
        if(tg_dNdpT_BW_ana_pid[i_mass]) delete tg_dNdpT_BW_ana_pid[i_mass];
        tg_dNdpT_BW_ana_pid[i_mass] = new TGraph();
        Double_t integral_BW = 0.0;
        for(Int_t i_pT = 0; i_pT < N_points_BW_ana; i_pT++)
        {
            Double_t pt_BW = i_pT*delta_pT + 0.0;
            Double_t v2_BW = 0;
            Double_t inv_yield_BW = 0;
            blastwave_yield_and_v2(pt_BW, arr_quark_mass_meson[i_mass], T_BW, Rho0_BW, Rho2_BW, RxOverRy_BW, inv_yield_BW, v2_BW);
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
            }
        }
        tg_dNdpT_BW_ana_pid[i_mass] -> SetLineColor(kAzure-2);
        tg_dNdpT_BW_ana_pid[i_mass] -> SetLineWidth(3);
        tg_dNdpT_BW_ana_pid[i_mass] -> SetLineStyle(1);
        if(!(fCheckBox_pid[i_mass]->GetState() == kButtonDown)) tg_dNdpT_BW_ana_pid[i_mass] -> SetLineStyle(9);
        if(fCheckBox_sel[1]->GetState() == kButtonDown)
        {
            tg_dNdpT_BW_ana_pid[i_mass] -> Draw("same L");
        }
    }

    printf("dNdpT ana plotted \n");
    fCanvasB->GetCanvas()->Modified();
    fCanvasB->GetCanvas()->Update();
    //------------------------------------------------------
}
//______________________________________________________________________________



//______________________________________________________________________________
void TBlastWaveGUI::MakePlotv2()
{
    printf("TBlastWaveGUI::MakePlotv2() \n");


    Int_t i_Temp   = vec_slider[0]->GetPosition();
    Int_t i_rho_0  = vec_slider[1]->GetPosition();
    Int_t i_rho_a  = vec_slider[2]->GetPosition();
    Int_t i_R_x    = vec_slider[3]->GetPosition();
    Int_t i_fboost = vec_slider[4]->GetPosition();

    //--------------------------------------------------------------------
    Double_t min_max_pT_range_pid_plot_ana[2][8];
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
    //--------------------------------------------------------------------

    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_make_plot_v2->ChangeBackground(green);


    //------------------------------------------------------------------------------------------------------------------------------------
    Double_t x_range_plot[2][2] =
    {
        {-0.2,4.3},
        {-0.2,14.7}
    };

    Double_t y_range_plot[2][2] =
    {
        {-0.12,0.34},
        {-0.12,0.34}
    };

    Int_t arr_color_plot[2][4] =
    {
        {kBlack,kGreen+1,kBlue+1,kRed+1},
        {kBlack,kGreen+1,kBlue+1,kRed+1}
    };
    Double_t arr_size_plot[2][4] =
    {
        {1.2,1.3,1.5,1.4},
        {1.2,1.3,1.5,1.4}
    };
    Int_t arr_marker_styleAB[2][4] =
    {
        {25,46,30,28},
        {25,46,30,28}
    };
    Int_t arr_marker_styleAB_x1[2][4] =
    {
        {21,47,29,34},
        {21,47,29,34}
    };
    Int_t arr_plot_orderAB[2][4]   =
    {
        {0,1,2,4},
        {3,5,6,7}
    };

    if(!c_2X1)
    {
        c_2X1 = new TCanvas("c_2X1","c_2X1",100,200,1500,620);
        c_2X1->SetTopMargin(0.02);
        c_2X1->SetBottomMargin(0.18);
        c_2X1->SetRightMargin(0.2);
        c_2X1->SetLeftMargin(0.2);

        c_2X1->SetLogy(0);
        c_2X1->Divide(2,1,0.0,0.0); // x divide, y divide, x margin, y margin
    }
    for(Int_t iPad = 0; iPad < 2; iPad++)
    {
        Double_t scaling_factor_2X1 = 0.78;
        Double_t Label_size_2X1     = 0.08;
        if(iPad == 0)
        {
            scaling_factor_2X1 = 0.78;
        }

        c_2X1->cd(iPad+1)->SetTicks(1,1);
        c_2X1->cd(iPad+1)->SetGrid(0,0);
        c_2X1->cd(iPad+1)->SetFillColor(10);
        c_2X1->cd(iPad+1)->SetRightMargin(0.01);
        c_2X1->cd(iPad+1)->SetTopMargin(0.01);
        HistName = "h_frame_2X1_";
        HistName += iPad;

        if(h_frame_2X1[iPad]) delete h_frame_2X1[iPad];
        h_frame_2X1[iPad] = c_2X1->cd(iPad+1)->DrawFrame(x_range_plot[iPad][0],y_range_plot[iPad][0],x_range_plot[iPad][1],y_range_plot[iPad][1],HistName.Data());
        h_frame_2X1[iPad]->SetStats(0);
        h_frame_2X1[iPad]->SetTitle("");
        h_frame_2X1[iPad]->GetXaxis()->SetTitleOffset(0.95/scaling_factor_2X1);
        h_frame_2X1[iPad]->GetYaxis()->SetTitleOffset(0.95/scaling_factor_2X1);
        h_frame_2X1[iPad]->GetXaxis()->SetLabelOffset(0.015*scaling_factor_2X1);
        //if(iPad == 0) h_frame_2X1[iPad]->GetXaxis()->SetLabelOffset(0.0*scaling_factor_2X1);
        h_frame_2X1[iPad]->GetYaxis()->SetLabelOffset(0.01*scaling_factor_2X1);
        h_frame_2X1[iPad]->GetXaxis()->SetLabelSize(Label_size_2X1*scaling_factor_2X1);
        h_frame_2X1[iPad]->GetYaxis()->SetLabelSize(Label_size_2X1*scaling_factor_2X1);
        h_frame_2X1[iPad]->GetXaxis()->SetTitleSize(Label_size_2X1*scaling_factor_2X1);
        h_frame_2X1[iPad]->GetYaxis()->SetTitleSize(Label_size_2X1*scaling_factor_2X1);
        h_frame_2X1[iPad]->GetXaxis()->SetNdivisions(505,'N');
        h_frame_2X1[iPad]->GetYaxis()->SetNdivisions(505,'N');
        h_frame_2X1[iPad]->GetXaxis()->CenterTitle();
        h_frame_2X1[iPad]->GetYaxis()->CenterTitle();
        h_frame_2X1[iPad]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        h_frame_2X1[iPad]->GetYaxis()->SetTitle("v_{2}");
        h_frame_2X1[iPad]->GetXaxis()->SetRangeUser(x_range_plot[iPad][0],x_range_plot[iPad][1]);
        h_frame_2X1[iPad]->GetYaxis()->SetRangeUser(y_range_plot[iPad][0],y_range_plot[iPad][1]);

        if(TL_line_base_plot[iPad]) delete TL_line_base_plot[iPad];
        TL_line_base_plot[iPad] = PlotLine(x_range_plot[iPad][0],x_range_plot[iPad][1],0.0,0.0,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)

        Int_t i_particle_use = 0;
        for(Int_t i_mass_loop = 0; i_mass_loop < 4; i_mass_loop++)
        {
            Int_t i_mass = arr_plot_orderAB[iPad][i_mass_loop];
            if(fCheckBox_pid[i_mass]->GetState() == kButtonDown)
            {
                //printf("Draw i_mass: %d \n", i_mass);
                tgae_v2_vs_pT_mesons_data_copy[i_mass] ->SetMarkerColor(kGray+1);
                tgae_v2_vs_pT_mesons_data_copy[i_mass] ->SetMarkerSize(arr_size_plot[iPad][i_mass_loop]*1.4);
                tgae_v2_vs_pT_mesons_data_copy[i_mass] ->SetMarkerStyle(arr_marker_styleAB_x1[iPad][i_mass_loop]);
                tgae_v2_vs_pT_mesons_data_copy[i_mass] ->Draw("same P");

                tgae_v2_vs_pT_mesons_data_copyB[i_mass] ->SetMarkerColor(kWhite);
                tgae_v2_vs_pT_mesons_data_copyB[i_mass] ->SetMarkerSize(arr_size_plot[iPad][i_mass_loop]*1.0);
                tgae_v2_vs_pT_mesons_data_copyB[i_mass] ->SetMarkerStyle(arr_marker_styleAB_x1[iPad][i_mass_loop]);
                tgae_v2_vs_pT_mesons_data_copyB[i_mass] ->Draw("same P");

                tgae_v2_vs_pT_mesons_data[i_mass] ->SetMarkerColor(arr_color_plot[iPad][i_mass_loop]);
                tgae_v2_vs_pT_mesons_data[i_mass] ->SetMarkerSize(arr_size_plot[iPad][i_mass_loop]);
                tgae_v2_vs_pT_mesons_data[i_mass] ->SetMarkerStyle(arr_marker_styleAB[iPad][i_mass_loop]);
                tgae_v2_vs_pT_mesons_data[i_mass] ->Draw("same P");


                tg_label_plot[0][i_mass] ->SetMarkerColor(kGray+1);
                tg_label_plot[0][i_mass] ->SetMarkerSize(arr_size_plot[iPad][i_mass_loop]*1.4);
                tg_label_plot[0][i_mass] ->SetMarkerStyle(arr_marker_styleAB_x1[iPad][i_mass_loop]);

                tg_label_plot[1][i_mass] ->SetMarkerColor(kWhite);
                tg_label_plot[1][i_mass] ->SetMarkerSize(arr_size_plot[iPad][i_mass_loop]*1.0);
                tg_label_plot[1][i_mass] ->SetMarkerStyle(arr_marker_styleAB_x1[iPad][i_mass_loop]);

                tg_label_plot[2][i_mass] ->SetMarkerColor(arr_color_plot[iPad][i_mass_loop]);
                tg_label_plot[2][i_mass] ->SetMarkerSize(arr_size_plot[iPad][i_mass_loop]);
                tg_label_plot[2][i_mass] ->SetMarkerStyle(arr_marker_styleAB[iPad][i_mass_loop]);

                Double_t x_pos_legend = 0.27 + iPad*0.57;
                Double_t y_pos_legend = 0.91  - i_particle_use*0.055;

                tg_label_plot[0][i_mass] ->SetPoint(0,0.1+iPad*11.9,0.3-i_particle_use*0.032);
                tg_label_plot[1][i_mass] ->SetPoint(0,0.1+iPad*11.9,0.3-i_particle_use*0.032);
                tg_label_plot[2][i_mass] ->SetPoint(0,0.1+iPad*11.9,0.3-i_particle_use*0.032);
                tg_label_plot[0][i_mass] ->Draw("same P");
                tg_label_plot[1][i_mass] ->Draw("same P");
                tg_label_plot[2][i_mass] ->Draw("same P");

                if(TLatex_legend_v2_plot[i_mass]) delete TLatex_legend_v2_plot[i_mass];
                TLatex_legend_v2_plot[i_mass] = plotTopLegend((char*)label_pid_spectra[i_mass].Data(),x_pos_legend,y_pos_legend,0.055,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
                i_particle_use++;

                if(fCheckBox_sel[0]->GetState() == kButtonDown)
                {
                    tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineWidth(3);
                    tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->SetLineColor(arr_color_plot[iPad][i_mass_loop]);
                    tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->Draw("same L hist");
                }
            }
        }


       
        //tp_v2_vs_pT_mesons[i_mass][i_R_x][i_fboost][i_Temp][i_rho_0][i_rho_a] ->GetXaxis()->SetRangeUser(min_val_pT,max_val_pT);

        const Int_t N_points_BW_ana = 35;
        printf("T_BW_fit_ana: %4.3f \n",T_BW_fit_ana);
        if(T_BW_fit_ana > 0.0)
        {
            for(int i_mass_loop = 0; i_mass_loop < 4; ++i_mass_loop)
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

                        blastwave_yield_and_v2(pt_BW, arr_quark_mass_meson[i_mass], T_BW_fit_ana, Rho0_BW_fit_ana, Rho2_BW_fit_ana, RxOverRy_BW_fit_ana, inv_yield_BW, v2_BW);
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
                    if(fCheckBox_sel[1]->GetState() == kButtonDown)
                    {
                        tg_v2_BW_ana_pid_plot[i_mass]       -> Draw("same L");
                        tg_v2_BW_ana_pid_plot_range[i_mass] -> Draw("same L");
                    }
                }
            }
        }
    }

    c_2X1 ->Modified();
    c_2X1 ->Update();
    printf("v2 ana plotted \n");
    //------------------------------------------------------------------------------------------------------------------------------------

}
//______________________________________________________________________________


//______________________________________________________________________________
void TBlastWaveGUI::MakePlotdNdpT()
{
    printf("TBlastWaveGUI::MakePlotdNdpT() \n");

    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_make_plot_dNdpT->ChangeBackground(green);
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
    if(c_2X1)
    {
        c_2X1 ->Write();
        c_2X1 ->SaveAs("v2_vs_pT_BW_fits.png");
    }

}
//______________________________________________________________________________


void v2_slider()
{
    new TBlastWaveGUI();
}
