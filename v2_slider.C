
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
    TGHorizontalFrame *hframeD2;
    //TGVerticalFrame *hframeD3;
    TGHorizontalFrame *hframeD3;
    TGVerticalFrame *hVframeD3;
    TGVerticalFrame *vframeD1;
    TGCompositeFrame *cframe2;
    TGTextButton *ButtonD1a;
    TGTextButton *Button_exit;
    TGTextButton *Button_save;
    TGHProgressBar* fHProg1 = NULL;
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

    TGHorizontalFrame* arr_HFrame_NEntry_limits[8];
    TGVerticalFrame* arr_VFrame_NEntry_limits[2];
    TGLabel*           arr_Label_NEntry_limits[2][8];
    TGNumberEntry*     arr_NEntry_limits[2][8];
    TGLayoutHints* arr_fL1[2];
    Double_t min_max_pT_range_pid[2][8];

    Double_t var_test = 5.2;

    TGTextButton *Button_minimize;
    TGTextButton *Button_stop_minimize;

    TFile* inputfile;
    TH1D* h_dummy;
    TH1D* h_dummy_dNdpT;
    TLegend* leg_v2_vs_pT_A = NULL;
    TLegend* leg_v2_vs_pT_B = NULL;
    Int_t flag_stop_minimize = 0;

    Double_t chi2_min;
    TString HistName;
    char NoP[50];

public:
    TTripleSliderDemo();
    virtual ~TTripleSliderDemo();
    void CloseWindow();
    void DoText(const char *text);
    void DoSlider();
    void DoMinimize();
    void StopMinimize();
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

    var_test = 2.6;
    //--------------------------------------------------------------------


    //--------------------------------------------------------------------
    min_max_pT_range_pid[0][0] = 0.1; // pi
    min_max_pT_range_pid[1][0] = 1.8;
    min_max_pT_range_pid[0][1] = 0.1; // K
    min_max_pT_range_pid[1][1] = 2.0;
    min_max_pT_range_pid[0][2] = 0.1; // p
    min_max_pT_range_pid[1][2] = 2.5;
    min_max_pT_range_pid[0][3] = 0.1; // phi
    min_max_pT_range_pid[1][3] = 2.5;
    min_max_pT_range_pid[0][4] = 0.1; // Omega
    min_max_pT_range_pid[1][4] = 2.8;
    min_max_pT_range_pid[0][5] = 0.1; // D0
    min_max_pT_range_pid[1][5] = 3.5;
    min_max_pT_range_pid[0][6] = 0.1; // J/Psi
    min_max_pT_range_pid[1][6] = 15.0;
    min_max_pT_range_pid[0][7] = 0.1; // Upsilon
    min_max_pT_range_pid[1][7] = 15.0;
    //--------------------------------------------------------------------


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
    //------------------------------------------------------------



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
    Button_save = new TGTextButton(hframeD1, "&Save ","gApplication->Terminate(0)");
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
    //fHProg1->Reset();
    cout << "A fHProg1: "<< fHProg1 << endl;
    //--------------



    //--------------
    // A horizontal frame
    hframeD2  = new TGHorizontalFrame(FrameD,200,100);

    // exit button
    Button_minimize = new TGTextButton(hframeD2, "Minimize ",10);
    Button_minimize->Connect("Clicked()", "TTripleSliderDemo", this, "DoMinimize()");
    hframeD2->AddFrame(Button_minimize, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // save button
    Button_stop_minimize = new TGTextButton(hframeD2, "Stop minimize ",10);
    Button_stop_minimize->Connect("Clicked()", "TTripleSliderDemo", this, "StopMinimize()");
    hframeD2->AddFrame(Button_stop_minimize, new TGLayoutHints(kLHintsCenterX,5,5,3,4));


    FrameD ->AddFrame(hframeD2, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------


    FrameD ->MapSubwindows();
    FrameD ->MapWindow();
    FrameD ->Resize(450,400); // size of frame
    FrameD ->Move(1250,750); // position of frame



    //--------------
    LHintsD4a = new TGLayoutHints(kLHintsCenterX,5,5,3,4);
    TGGC myGC = *gClient->GetResourcePool()->GetFrameGC();
    TGFont *myfont = gClient->GetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");

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
    cout << "DoSlider started" << endl;

    Double_t number_entry = arr_NEntry_limits[0][0]->GetNumberEntry()->GetNumber();
    printf("number_entry: %4.1f \n",number_entry);

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

        HistName = "#chi^{2}/ndf = ";
        sprintf(NoP,"%4.2f",chi2_min);
        HistName += NoP;
        plotTopLegend((char*)HistName.Data(),0.18,0.85,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

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
    
    Double_t chi2[3]; 
    
    for(Int_t i_mass = 0; i_mass < 3; i_mass++)
    {
        chi2[i_mass] = 0.0;

        if(tp_v2_vs_pT_mesons[i_mass])
        {
            Int_t plot_centrality   = 4;
            Int_t n_arr             = vec_graphs[plot_centrality+14*i_mass]->GetN();
            Double_t x_pid;
            Double_t v2_pid;
            Double_t v2_bw_pid;

            for(Int_t i_pT = 0; i_pT < n_arr; i_pT++) // pT loop
            {
                // Loop over pT points
                //tp_v2_vs_pT_mesons[i_mass] ->GetPoint...
                // Compare to data vec_graphs[plot_centrality]...

                vec_graphs[plot_centrality+14*i_mass]            ->GetPoint(i_pT,x_pid,v2_pid);
                v2_bw_pid           = tp_v2_vs_pT_mesons[i_mass] ->GetBinContent(tp_v2_vs_pT_mesons[i_mass]->FindBin(x_pid));
                chi2[i_mass]        += (v2_pid-v2_bw_pid)*(v2_pid-v2_bw_pid);

                //if(i_pid == 0) printf("i_pid: %d, i_cent: %d, i_pT: %d, pT: %4.2f, v2: %4.3f, dNdpT: %4.6f, v2_unnorm_pid: %4.3f \n",i_pid,i_cent,i_pT,x_arr_pid[i_pT],y_arr_pid,y_pt_arr_pid,v2_unnorm_pid[i_pT]);
            }

            cout << "i_mass: " << i_mass << endl;
            cout << "chi2[i_mass]: " << chi2[i_mass] << endl;
        }
    }


    //if(fHProg1) cout << "OK" << endl;
    //else cout << "Not OK" << endl;

}

//______________________________________________________________________________
void TTripleSliderDemo::StopMinimize()
{
    printf("Minimization stopped \n");
    flag_stop_minimize = 1;

}

//______________________________________________________________________________
void TTripleSliderDemo::DoMinimize()
{
    cout << "DoMinimize started" << endl;

    flag_stop_minimize = 0;

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

    chi2_min = 10000000.0;
    Double_t chi2_tot = 0;

    Double_t N_total_params = TMath::Power(8,5);

    Double_t N_params_use = 0;
    Int_t    N_params_total_use = 0;
    Double_t fraction_progress_bar_update = 0.005;

    Int_t n_arr;
    Int_t arr_mass[8] = {0,1,2,6,7}; // mapping array, need to be changed later on

    for (Int_t i_Temp = 0; i_Temp < 8; i_Temp++)
    {
        if(flag_stop_minimize) break;
        //printf("i_Temp: %d \n",i_Temp);
        for (Int_t i_rho_0 = 0; i_rho_0 < 8; i_rho_0++)
        {
            if(flag_stop_minimize) break;
            //printf("i_rho_0: %d \n",i_rho_0);
            for (Int_t i_rho_a = 0; i_rho_a < 8; i_rho_a++)
            {
                if(flag_stop_minimize) break;
                //printf("i_rho_a: %d \n",i_rho_a);
                for (Int_t i_R_x = 0; i_R_x < 8; i_R_x++)
                {
                    if(flag_stop_minimize) break;
                    for (Int_t i_fboost = 0; i_fboost < 8; i_fboost++)
                    {
                        if(flag_stop_minimize) break;

                        //chi2 - sum - the one we need
                        //Double_t chi2[5]; //chi2 - individual - in case we want to look at it at some point
                        Int_t    nop_tot = 0; // number of points - sum
                        Double_t pT_lim[5] = {1.8,2,2.5,TMath::MaxElement(tg_JPsi_v2_vs_pT->GetN(),tg_JPsi_v2_vs_pT->GetX()),TMath::MaxElement(tg_Upsilon_v2_vs_pT->GetN(),tg_Upsilon_v2_vs_pT->GetX())};


                        

                        for(Int_t i_mass = 0; i_mass < 5; i_mass++)
                        {
                            Int_t i_mass_arr = arr_mass[i_mass];

                            tp_v2_vs_pT_mesons[i_mass] = NULL;

                            tp_v2_vs_pT_mesons[i_mass] = (TProfile*)inputfile->Get(Form("v2_vs_pT_BW_id%d_T%d_rho0%d_rhoa%d_Rx%d_fb%d",i_mass,i_Temp,i_rho_0,i_rho_a,i_R_x,i_fboost));


                            //chi2[i_mass] = 0.0; // for everybody
                            n_arr        = 0;

                            if(tp_v2_vs_pT_mesons[i_mass])
                            {
                                Int_t plot_centrality   = 4;

                                // get n_arr for different particles

                                if(i_mass < 3)
                                {
                                    n_arr             = vec_graphs[plot_centrality+14*i_mass]->GetN();
                                }

                                if(i_mass == 3)
                                {
                                    n_arr             = tg_JPsi_v2_vs_pT                     ->GetN();
                                }

                                if (i_mass == 4)
                                {
                                    n_arr             = tg_Upsilon_v2_vs_pT                  ->GetN();
                                }

                                Double_t x_pid; //for everybody too
                                Double_t v2_pid;
                                Double_t v2_err_pid;
                                Double_t v2_bw_pid;
                                //Int_t nop = 0; // number of points - individual

                                for(Int_t i_pT = 0; i_pT < n_arr; i_pT++) // pT loop
                                {
                                    //get v2_pid, v2_err_pid for different particles
                                    if(i_mass < 3)
                                    {
                                        vec_graphs[plot_centrality+14*i_mass]                       ->GetPoint(i_pT,x_pid,v2_pid);
                                        v2_err_pid          = vec_graphs[plot_centrality+14*i_mass] ->GetErrorY(i_pT);
                                    }

                                    if(i_mass == 3)
                                    {
                                        tg_JPsi_v2_vs_pT                       ->GetPoint(i_pT,x_pid,v2_pid);
                                        v2_err_pid          = tg_JPsi_v2_vs_pT ->GetErrorY(i_pT);
                                    }

                                    if(i_mass == 4)
                                    {
                                        tg_Upsilon_v2_vs_pT                       ->GetPoint(i_pT,x_pid,v2_pid);
                                        v2_err_pid          = tg_Upsilon_v2_vs_pT ->GetErrorY(i_pT);
                                    }
                                        
                                    //if(x_pid <= pT_lim[i_mass] && v2_err_pid != 0) // calculate only within cetrain pT range; one loop for everybody
                                    if(x_pid >= min_max_pT_range_pid[0][i_mass_arr]
                                       && x_pid <= min_max_pT_range_pid[1][i_mass_arr]
                                       && v2_err_pid != 0) // calculate only within cetrain pT range; one loop for everybody
                                    {
                                        v2_bw_pid           = tp_v2_vs_pT_mesons[i_mass]        ->GetBinContent(tp_v2_vs_pT_mesons[i_mass]->FindBin(x_pid));
                                        //chi2[i_mass]        += ((v2_pid-v2_bw_pid)*(v2_pid-v2_bw_pid))/(v2_err_pid*v2_err_pid);
                                        chi2_tot            += ((v2_pid-v2_bw_pid)*(v2_pid-v2_bw_pid))/(v2_err_pid*v2_err_pid);
                                        //nop                 += 1;
                                        nop_tot             += 1;
                                     }

                                     //cout << "i_pT: " << i_pT << endl;

                                }

                                
                                // if (nop > 5) //individual chi2 
                                // {
                                //     chi2[i_mass] = chi2[i_mass]/(nop-5);
                                // }

                            }

                        }


                        if (nop_tot > 5) //sum of chi2
                        {
                            chi2_tot = chi2_tot/(nop_tot-5);
                        }

                        //chi2_min = chi_tot;
                        //cout << "vec_graphs[i_cent]->GetN():" << vec_graphs[i_cent+i_pid*7+pid_helper]->GetN() << endl;
                        //cout << "i_pid:" << i_pid << endl;
                        if(chi2_tot < chi2_min)
                        {
                            chi2_min     = chi2_tot;
                            i_Temp_min   = i_Temp;
                            i_rho_0_min  = i_rho_0;
                            i_rho_a_min  = i_rho_a;
                            i_R_x_min    = i_R_x;
                            i_fboost_min = i_fboost;

                            printf("chi2_min: %4.3f \n",chi2_min);

                            vec_slider[0]->SetPosition(i_Temp_min);
                            vec_slider[1]->SetPosition(i_rho_0_min);
                            vec_slider[2]->SetPosition(i_rho_a_min);
                            vec_slider[3]->SetPosition(i_R_x_min);
                            vec_slider[4]->SetPosition(i_fboost_min);
                            DoSlider();
                            gSystem->Sleep(100);
                            gSystem->ProcessEvents();
                        }

                        //cout << "chi2_min: " << chi2_m << endl;
                        //cout << "chi2_tot: " << chi2_tot << endl;



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

    Double_t Temp_loop_start  = 0.08;
    Double_t rho_0_loop_start = 0.3;
    Double_t arr_rho_a[8]   = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35};
    Double_t arr_R_x[8]     = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
    Double_t arr_f_boost[8] = {0.05,0.1,0.15,0.2,0.4,0.6,0.8,1.0};

    Double_t Temp_val_min   = Temp_loop_start  + i_Temp_min*0.02;
    Double_t rho_0_val_min  = rho_0_loop_start + i_rho_0_min*0.125;
    Double_t rho_a_val_min  = arr_rho_a[i_rho_a_min];
    Double_t R_x_val_min    = arr_R_x[i_R_x_min];
    Double_t fboost_val_min = arr_f_boost[i_fboost_min];

    Double_t arr_param_val_min[5] = {Temp_val_min,rho_0_val_min,rho_a_val_min,R_x_val_min,fboost_val_min};
    

    cout << "Hello Dave. Your minimum Chi2 = " << chi2_min << endl;
    cout << "Corresponding parameters are:" << chi2_min << endl;
    cout << "Temperature: " << arr_param_val_min[0] << endl;
    cout << "rho_0: " << arr_param_val_min[1] << endl;
    cout << "rho_a: " << arr_param_val_min[2] << endl;
    cout << "R_x: " << arr_param_val_min[3] << endl;
    cout << "fboost: " << arr_param_val_min[4] << endl;

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
