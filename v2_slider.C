
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
class TTripleSliderDemo : public TGMainFrame {
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
    TGCheckButton       *fCheck1, *fCheck2, *fCheck3;

    vector<TProfile*> tp_v2_vs_pT_mesons;
    vector<TH1D*>     h_dN_dpT_mesons;
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
    TGVerticalFrame *vframeD1;
    TGCompositeFrame *cframe2;
    TGTextButton *ButtonD1a;
    TGTextButton *Button_exit;
    TGTextButton *Button_save;

    TFile* inputfile;
    TH1D* h_dummy;
    TH1D* h_dummy_dNdpT;
    TLegend* leg_v2_vs_pT_A = NULL;
    TLegend* leg_v2_vs_pT_B = NULL;

public:
    TTripleSliderDemo();
    virtual ~TTripleSliderDemo();
    void CloseWindow();
    void DoText(const char *text);
    void DoSlider();
    void HandleButtons();
    ClassDef(TTripleSliderDemo, 0)
};

//______________________________________________________________________________
TTripleSliderDemo::TTripleSliderDemo() : TGMainFrame(gClient->GetRoot(), 100, 100)
{
    //--------------------------------------------------------------------
    cout << "v2_slider started" << endl;

    TGaxis::SetMaxDigits(3);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //--------------------------------------------------------------------

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

    inputfile = TFile::Open("./Data/merge_v2_boost.root");
    tp_v2_vs_pT_mesons.resize(5);
    h_dN_dpT_mesons.resize(5);

    init_data();
    init_pT_spectra_data();
    init_JPsi_spectra_data();
    make_5_60_spectra();
    Init_v2_Mathematica();

    char buf[32];
    SetCleanup(kDeepCleanup);

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
        vec_slider[i_param]->Connect("PositionChanged(Int_t)", "TTripleSliderDemo",this, "DoSlider()");
        vec_slider[i_param]->SetRange(0,7);
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
        vec_TextEntry[i_param]->Connect("TextChanged(char*)", "TTripleSliderDemo", this,"DoText(char*)");
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


    FrameB = new TGMainFrame(gClient->GetRoot(), 400, 100);
    FrameB ->SetWindowName("Transverse momentum spectra");
    FrameB ->MapSubwindows();
    FrameB ->Resize(GetDefaultSize());
    FrameB ->MapWindow();
    FrameB->Resize(1200,700);


    // Create an embedded canvas and add to the main frame, centered in x and y
    // and with 30 pixel margins all around
    fCanvasB = new TRootEmbeddedCanvas("CanvasB", FrameB, 1200, 700);
    fLcanC = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 10);
    FrameB ->AddFrame(fCanvasB, fLcanC);
    fCanvasB->GetCanvas()->Divide(3,2);
    for(Int_t iPad = 1; iPad <= 6; iPad++)
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
    FrameB ->Move(700,10);

    for(Int_t iPad = 1; iPad <= 6; iPad++)
    {
        fCanvasB ->GetCanvas()->cd(iPad);
        h_dummy  ->Draw();
        fCanvasB ->GetCanvas()->cd(iPad)->Modified();
        fCanvasB ->GetCanvas()->cd(iPad)->Update();
    }


    FrameD = new TGMainFrame(gClient->GetRoot(), 400, 100);
    FrameD ->SetWindowName("Buttons");
    FrameD ->MapSubwindows();
    FrameD ->Resize(GetDefaultSize());
    FrameD ->MapWindow();
    FrameD ->Resize(400,400); // size of frame

    //hframeD1  = new TGHorizontalFrame(FrameD,200,100);
    vframeD1 = new TGVerticalFrame(FrameD, 150, 150);
    cframe2 = new TGCompositeFrame(vframeD1, 170, 50,kHorizontalFrame | kFixedWidth);

    // exit button
    Button_exit = new TGTextButton(cframe2, "&Exit ","gApplication->Terminate(0)");
    cframe2->AddFrame(Button_exit, new TGLayoutHints(kLHintsTop | kLHintsExpandX,3, 2, 2, 2));

    // save button
    Button_save = new TGTextButton(cframe2, "&Save ","gApplication->Terminate(0)");
    cframe2->AddFrame(Button_save, new TGLayoutHints(kLHintsTop | kLHintsExpandX,2, 0, 2, 2));


    vframeD1->AddFrame(cframe2, new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 1));

    // A button
    //ButtonD1a = new TGTextButton(hframeD1,"&DrawD1a");
    //hframeD1->AddFrame(ButtonD1a, new TGLayoutHints(kLHintsCenterX,5,5,3,4)); // left, right, top, bottom

    FrameD->AddFrame(vframeD1, new TGLayoutHints(kLHintsExpandX, 2, 2, 5, 1));

    // Set a name to the main frame
    FrameD->SetWindowName("A frame with buttons");

    // Map all subwindows of main frame
    FrameD->MapSubwindows();

    // Map main frame
    FrameD->MapWindow();
    FrameD ->Move(1250,750); // position of frame



    DoSlider();

}
//______________________________________________________________________________
TTripleSliderDemo::~TTripleSliderDemo()
{
    // Clean up
    Cleanup();
}
//______________________________________________________________________________
void TTripleSliderDemo::CloseWindow()
{
    // Called when window is closed via the window manager.
    delete this;
}


//______________________________________________________________________________
void TTripleSliderDemo::DoText(const char * /*text*/)
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
void TTripleSliderDemo::DoSlider()
{

    //0.18, 1.175, 0.15, 0.7, 0.05
    // 0.18, 0.8, 0.15, 0.7, 0.05

    //cout << "DoSlider()" << endl;

    Int_t i_Temp   = vec_slider[0]->GetPosition();
    Int_t i_rho_0  = vec_slider[1]->GetPosition();
    Int_t i_rho_a  = vec_slider[2]->GetPosition();
    Int_t i_R_x    = vec_slider[3]->GetPosition();
    Int_t i_fboost = vec_slider[4]->GetPosition();

    Double_t Temp_loop_start  = 0.08;
    Double_t rho_0_loop_start = 0.3;
    Double_t arr_rho_a[8]   = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35};
    Double_t arr_R_x[8]     = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
    Double_t arr_f_boost[8] = {0.05,0.1,0.15,0.2,0.4,0.6,0.8,1.0};

    Double_t Temp_val   = Temp_loop_start  + i_Temp*0.02;
    Double_t rho_0_val  = rho_0_loop_start + i_rho_0*0.125;
    Double_t rho_a_val  = arr_rho_a[i_rho_a];
    Double_t R_x_val    = arr_R_x[i_R_x];
    Double_t fboost_val = arr_f_boost[i_fboost];

    Double_t arr_param_val[5] = {Temp_val,rho_0_val,rho_a_val,R_x_val,fboost_val};

    Int_t arr_color_line[5] = {kBlack,kGreen,kRed,kMagenta,kCyan};

    for(Int_t i_mass = 0; i_mass < 5; i_mass++)
    {
        tp_v2_vs_pT_mesons[i_mass] = NULL;
        tp_v2_vs_pT_mesons[i_mass] = (TProfile*)inputfile->Get(Form("v2_vs_pT_BW_id%d_T%d_rho0%d_rhoa%d_Rx%d_fb%d",i_mass,i_Temp,i_rho_0,i_rho_a,i_R_x,i_fboost));

        if(tp_v2_vs_pT_mesons[i_mass])
        {
            tp_v2_vs_pT_mesons[i_mass] ->SetLineWidth(3);
            tp_v2_vs_pT_mesons[i_mass] ->SetLineColor(arr_color_line[i_mass]);
            tp_v2_vs_pT_mesons[i_mass] ->SetLineStyle(1);
        }

        h_dN_dpT_mesons[i_mass] = NULL;
        h_dN_dpT_mesons[i_mass] = (TH1D*)inputfile->Get(Form("h_dN_dpT_vs_pT_BW_id%d_T%d_rho0%d_rhoa%d_Rx%d_fb%d",i_mass,i_Temp,i_rho_0,i_rho_a,i_R_x,i_fboost));

        if(h_dN_dpT_mesons[i_mass])
        {
            h_dN_dpT_mesons[i_mass] ->SetLineWidth(3);
            h_dN_dpT_mesons[i_mass] ->SetLineColor(arr_color_line[i_mass]);
            h_dN_dpT_mesons[i_mass] ->SetLineStyle(1);
        }
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

    //fCanvasB ->GetCanvas() ->cd(1);
    //h_dummy ->Draw();


    //-------------------------------------
    fCanvas ->GetCanvas() ->cd();
    if(tp_v2_vs_pT_mesons[0])
    {
        h_dummy ->Draw();

        PlotLine(0.0,15.0,0.0,0.0,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)

        Int_t plot_centrality = 4;


        // 30-40%
        vec_graphs[plot_centrality] ->SetMarkerColor(arr_color_mass[0]);
        vec_graphs[plot_centrality] ->Draw("same P"); // pions
        //vec_graphs[plot_centrality+7] ->Draw("same P"); // charged kaons
        vec_graphs[plot_centrality+14] ->SetMarkerColor(arr_color_mass[1]);
        vec_graphs[plot_centrality+14] ->Draw("same P"); // K0s
        vec_graphs[plot_centrality+28] ->SetMarkerColor(arr_color_mass[2]);
        vec_graphs[plot_centrality+28] ->Draw("same P"); // protons


#if 0
    for(Int_t i_mass = 0; i_mass < 5; i_mass++)
    {
        vec_tg_v2_vs_pT_Mathematica[i_mass] ->SetMarkerColor(arr_color_mass[i_mass]);
        vec_tg_v2_vs_pT_Mathematica[i_mass] ->SetMarkerStyle(28); // 24
        vec_tg_v2_vs_pT_Mathematica[i_mass] ->SetLineColor(arr_color_mass[i_mass]);
        vec_tg_v2_vs_pT_Mathematica[i_mass] ->SetLineStyle(1);
        vec_tg_v2_vs_pT_Mathematica[i_mass] ->SetLineWidth(2);
        vec_tg_v2_vs_pT_Mathematica[i_mass] ->Draw("same P");
    }
#endif

        /*
        vec_tge_v2_vs_pT_560_pid[0] ->SetMarkerColor(arr_color_mass[0]);
        vec_tge_v2_vs_pT_560_pid[0] ->Draw("same P"); // pions
        vec_tge_v2_vs_pT_560_pid[1] ->SetMarkerColor(arr_color_mass[1]);
        vec_tge_v2_vs_pT_560_pid[1] ->Draw("same P"); // K0s
        vec_tge_v2_vs_pT_560_pid[2] ->SetMarkerColor(arr_color_mass[2]);
        vec_tge_v2_vs_pT_560_pid[2] ->Draw("same P"); // protons
        */

        tg_JPsi_v2_vs_pT    ->SetMarkerColor(arr_color_mass[3]);
        tg_JPsi_v2_vs_pT    ->Draw("same P");
        tg_D0_v2_vs_pT      ->SetMarkerColor(kGray+1);
        tg_D0_v2_vs_pT      ->Draw("same P");
        tg_Upsilon_v2_vs_pT ->SetMarkerColor(arr_color_mass[4]);
        tg_Upsilon_v2_vs_pT ->Draw("same P");

        tp_v2_vs_pT_mesons[0] ->GetXaxis()->SetRangeUser(0.0,4.0);
        tp_v2_vs_pT_mesons[1] ->GetXaxis()->SetRangeUser(0.0,4.0);
        tp_v2_vs_pT_mesons[2] ->GetXaxis()->SetRangeUser(0.0,4.0);
        tp_v2_vs_pT_mesons[0] ->SetLineWidth(4);
        tp_v2_vs_pT_mesons[1] ->SetLineWidth(4);
        tp_v2_vs_pT_mesons[2] ->SetLineWidth(4);
        tp_v2_vs_pT_mesons[3] ->SetLineWidth(4);
        tp_v2_vs_pT_mesons[4] ->SetLineWidth(4);
        tp_v2_vs_pT_mesons[0] ->Draw("same L hist");
        tp_v2_vs_pT_mesons[1] ->Draw("same L hist");
        tp_v2_vs_pT_mesons[2] ->Draw("same L hist");
        tp_v2_vs_pT_mesons[3] ->Draw("same L hist");
        tp_v2_vs_pT_mesons[4] ->Draw("same L hist");

        if(leg_v2_vs_pT_A) delete leg_v2_vs_pT_A;
        leg_v2_vs_pT_A = new TLegend(0.76,0.63,0.86,0.82); // x1,y1,x2,y2
        leg_v2_vs_pT_A->SetBorderSize(0);
        leg_v2_vs_pT_A->SetFillColor(0);
        leg_v2_vs_pT_A->SetTextSize(0.045);
        //leg_v2_vs_pT_A->AddEntry((TH1D*)vec_tge_v2_vs_pT_560_pid[0]->Clone(),"#pi","p");
        //leg_v2_vs_pT_A->AddEntry((TH1D*)vec_tge_v2_vs_pT_560_pid[1]->Clone(),"K_{s}^{0}","p");
        //leg_v2_vs_pT_A->AddEntry((TH1D*)vec_tge_v2_vs_pT_560_pid[2]->Clone(),"p","p");
        leg_v2_vs_pT_A->AddEntry((TH1D*)vec_graphs[plot_centrality]->Clone(),"#pi","p");
        leg_v2_vs_pT_A->AddEntry((TH1D*)vec_graphs[plot_centrality+14]->Clone(),"K_{s}^{0}","p");
        leg_v2_vs_pT_A->AddEntry((TH1D*)vec_graphs[plot_centrality+28]->Clone(),"p","p");
        leg_v2_vs_pT_A->Draw();

        if(leg_v2_vs_pT_B) delete leg_v2_vs_pT_B;
        leg_v2_vs_pT_B = new TLegend(0.87,0.68,0.97,0.81); // x1,y1,x2,y2
        leg_v2_vs_pT_B->SetBorderSize(0);
        leg_v2_vs_pT_B->SetFillColor(0);
        leg_v2_vs_pT_B->SetTextSize(0.045);
        leg_v2_vs_pT_B->AddEntry((TGraphErrors*)tg_Upsilon_v2_vs_pT->Clone(),"#Upsilon","p");
        leg_v2_vs_pT_B->AddEntry((TGraphErrors*)tg_JPsi_v2_vs_pT->Clone(),"J/#Psi","p");
        leg_v2_vs_pT_B->Draw();

        plotTopLegend((char*)"30-40%",0.73,0.83,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        plotTopLegend((char*)"5-60%",0.85,0.83,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        plotTopLegend((char*)"|y|<0.5",0.73,0.89,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        plotTopLegend((char*)"2.5<y<4",0.85,0.89,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        fCanvas->GetCanvas()->Modified();
        fCanvas->GetCanvas()->Update();
    }
    //-------------------------------------


    //-------------------------------------
    // pt spectra
    for(Int_t iPid = 0; iPid < 5; iPid++)
    {
        if(h_dN_dpT_mesons[iPid])
        {
            Double_t integral = h_dN_dpT_mesons[iPid] ->Integral("width");
            if(integral <= 0.0) continue;
            h_dN_dpT_mesons[iPid] ->Scale(1.0/integral);
        }
    }

    for(Int_t i_pid = 0; i_pid < 3; i_pid++)
    {
        vec_tgae_pT_spectra[i_pid][3]->SetMarkerColor(kBlue);
        vec_tgae_pT_spectra[i_pid][10]->SetMarkerColor(kGreen);
    }

    Double_t x_range_dNdpT[5] = {3.5,4.1,4.5,8.5,15.0};
    Double_t y_range_dNdpT[5] = {2.2,1.4,1.1,0.45,0.35};

    for(Int_t iPid = 0; iPid < 5; iPid++)
    {
        fCanvasB ->GetCanvas()->cd(iPid + 1);
        h_dummy_dNdpT ->GetXaxis()->SetRangeUser(0.0,x_range_dNdpT[iPid]);
        h_dummy_dNdpT ->GetYaxis()->SetRangeUser(-0.2,y_range_dNdpT[iPid]);
        h_dummy_dNdpT ->DrawCopy("");
        if(iPid < 3)
        {
            //vec_tgae_pT_spectra[iPid][10] ->Draw("same P"); // data, 5-60%
            vec_tgae_pT_spectra[iPid][7] ->Draw("same P"); // data, 30-40%
        }
        if(iPid == 3)
        {
            //tge_JPsi_spectra[1][0] ->SetLineColor(kRed);
            //tge_JPsi_spectra[1][0] ->Draw("same P");; // 0-20%, 20-40%, 40-90%
            tge_JPsi_forward_spectrum_stat ->SetLineColor(kGreen);
            tge_JPsi_forward_spectrum_stat ->SetMarkerColor(kGreen);
            tge_JPsi_forward_spectrum_stat ->SetMarkerStyle(20);
            tge_JPsi_forward_spectrum_stat ->Draw("same P");
        }
        if(iPid == 4)
        {
            // No Upsilon pT spectrum available
        }

        h_dN_dpT_mesons[iPid] ->SetLineColor(kRed); // blast wave MC
        h_dN_dpT_mesons[iPid] ->SetLineWidth(5); // blast wave MC
        h_dN_dpT_mesons[iPid] ->DrawCopy("same hist L"); // blast wave MC
    }

    for(Int_t iPad = 1; iPad <= 3; iPad++)
    {
        fCanvasB ->GetCanvas()->cd(iPad);
        plotTopLegend((char*)label_pid_spectra[iPad-1].Data(),0.75,0.83,0.06,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        plotTopLegend((char*)"|y|<0.5",0.75,0.77,0.06,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    }
    for(Int_t iPad = 4; iPad <= 5; iPad++)
    {
        fCanvasB ->GetCanvas()->cd(iPad);
        plotTopLegend((char*)label_pid_spectra[iPad-1].Data(),0.75,0.83,0.06,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        plotTopLegend((char*)"2.5<y<4",0.75,0.77,0.06,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    }

    fCanvasB->GetCanvas()->Modified();
    fCanvasB->GetCanvas()->Update();
    //-------------------------------------


    // Do chi2 calculation here:
    /*
    Double_t chi2 = 0.0;
    for(Int_t i_mass = 0; i_mass < 5; i_mass++)
    {
        if(tp_v2_vs_pT_mesons[i_mass])
        {
            // Loop over pT points
            //tp_v2_vs_pT_mesons[i_mass] ->GetPoint...
            // Compare to data vec_graphs[plot_centrality]...
            chi2 += ...
        }
    }
    */
}
//______________________________________________________________________________
void TTripleSliderDemo::HandleButtons()
{
    // Handle different buttons.
    TGButton *btn = (TGButton *) gTQSender;
    Int_t id = btn->WidgetId();
    switch (id) {
    case HCId1:
        fHslider1->SetConstrained(fCheck1->GetState());
        break;
    case HCId2:
        fHslider1->SetRelative(fCheck2->GetState());
        break;
    default:
        break;
    }
}

void v2_slider()
{
    new TTripleSliderDemo();
}
