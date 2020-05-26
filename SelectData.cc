
// Store hep data file in your local file
// Files button opens local file. Select root file & load it with Open button
// Choose output file name
// Add as much entries as needed with the Add entry button
// Choose between v2 or dNdpt, select particle, energy, centrality and the table with the data
// Click Plot button for plotting or Save button for saving the graphs in the output file


// Nadine Gruenwald 03/25/2020

// Add option to open root or text files. When openig text files you can switch between pt or pt_low and pt_high.
// Before saving you can select if dNdpt is divided by pt or not.
// 03/27/2020

// Can now switch between stat_error and syst_error or select both
// 03/30/2020


// Add combo box to select graph number (hep files)
// 04/17/2020

// Add function Multiply_pT. Multiply dNdpt with pt if checkbox dNdpt is clicked
// 04/20/2020

// Add check boxes mt. dNdmt, 1/mt*dNdmt
// 04/24/2020


// Add Lambda_c
// 05/25/2020

// Missing : plot dNdpt, dNdmT only for text files
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------

bool is_whitespace(const std::string& s)
{
  for (std::string::const_iterator it = s.begin(); it != s.end(); ++it) {
    if (!isspace(*it)) {
      return false;
    }
  }
  return true;
}

vector<Double_t> extract_values_from_string(std::string str, Int_t N_columns)
{
    std::size_t found[3];
    std::string sub_str;
    vector<Double_t> values;

    for(Int_t i_column = 0; i_column < (N_columns-1); i_column++)
    {
        found[0]  = str.find("\t");
        found[1]  = str.find(" ");
        if(found[1] < found[0]) found[0] = found[1];
        sub_str = str.substr(0,found[0]);
        //cout << "sub_str column: -" << sub_str << "-" << endl;
        values.push_back(std::stod(sub_str));
        str.replace(0,found[0]+1,"");

        //printf("i_column: %d, value: %4.6f \n",i_column,std::stod(sub_str));

        sub_str = "\t";
        bool is_empty = 1;
        Int_t counter = 0;
        while(sub_str == "\t" || sub_str == " " || is_empty)
        {
            found[0] = str.find_first_not_of("\t");
            found[1] = str.find_first_not_of(' ');
            if(found[1] > found[0]) found[0] = found[1];
            sub_str = str.substr(0,found[0]);
            str.replace(0,found[0],"");

            is_empty = is_whitespace(sub_str);
            //if(is_empty) cout << "Only empty" << endl;

            //cout << "counter: " << counter << ", sub_str: -" << sub_str << "-" << endl;
            //cout << "str: " << str << endl;

            counter++;
            if(counter > 20) break;
        }
    }
    values.push_back(std::stod(str));

    return values;
}


//Plot title of x- and y-axis
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1)
{

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

//Plot zero line
void PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
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
}



static TFile *newfile; // Input file
static TGRadioButton *hep, *text,*pt, *pt_high_low, *mt, *mt_low_high;
static std::fstream myfile;
Int_t flag_file_opened = 0;
static TGCheckButton *stat_error, *syst_error, *syst_percent, *syst_low_high, *pt_error;
//--------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------
class MyMainFrame : public TGMainFrame
{

private:
    TGMainFrame         *fMain;
    TGButton            *Exit, *Add_entry, *Remove_entry, *Add_file, *Save, *Plot;
    TGCheckButton       *dNdptpt, *dNdmt, *dNdmtmt;
    TGTextEntry         *Add_output_file_name;
    TGCompositeFrame    *hframebutton, *hframefiles, *hframeoutput,  *hframe, *V2frame, *PIDframe, *Energyframe, *Centralityframe, *Tableframe,*Graphframe, *Centrality_labelframe;
    TGRadioButton       *v2, *dNdpt;
    TGLayoutHints       *fL_main, *fL_framebutton, *fL_framefiles, *fL_button,  *fL_hframe, *fL_Energyframe, *fL_Centralityframe, *fL_upperlower;
    TGComboBox          *fCombo, *ComboPID;
    TGNumberEntry       *table, *graph, *energy_entry, *centrality_lower, *centrality_upper;
    TGVButtonGroup      *fButtonGroup;
    TGLabel             *label_PID, *label_Energy, *label_Table, *label_Centrality, *label_upper, *label_lower, *label_space;
    TCanvas             *c_3X3 = NULL;
    TFile               *RootFile = NULL;
    TH1D                *h_dummy;
    Int_t                N_used_entries = 0;
    TString              HistName;
    Double_t mass_pid[23] = {0.13957,0.13957,0.493677,0.493677,0.938272,0.938272,1.019460,1.32171, 1.32171, 1.67245,1.67245,1.115683,1.115683,2.28646,0.497611,1.86962,3.096916,9.46030,1.875612,1.875612, 2.8094313, 2.8094313, 2.80945};
    TString label_full_pid_spectra[23] = {"Pi+","Pi-","K+","K-","P","Pbar", "Phi","Xi-","Xibar+","Omega-","Omegabar+","Lambda","Lambdabar","Lambda_C", "K0S","D0", "J/Psi","Upsilon","d","dbar","He3", "He3bar", "t"}; // 9 -> 21




    vector<TGCompositeFrame*> vec_hframe, vec_V2frame, vec_PIDframe, vec_Tableframe, vec_Graphframe, vec_Centrality_labelframe, vec_Centralityframe, vec_Energyframe;
    vector<TGRadioButton*>    vec_v2, vec_dNdpt;
    vector<TGComboBox*>       vec_fCombo;
    vector<TGNumberEntry*>    vec_table, vec_graph, vec_energy_entry, vec_centrality_lower, vec_centrality_upper;
    vector<TGVButtonGroup*>   vec_fButtonGroup;
    vector<TGLabel*>          vec_label_PID, vec_label_Energy, vec_label_Table,vec_label_Graph, vec_label_Centrality, vec_label_upper, vec_label_lower, vec_label_space;
    vector<TString>           vec_HistName;
    vector<TGraphAsymmErrors*> vec_tgae ,vec_tgae_syst;

public:
    MyMainFrame(const TGWindow *p, UInt_t w, UInt_t );
    virtual ~MyMainFrame();

    //slots
    void CloseWindow();
    void v2_clicked(Int_t index_frame);
    void dNdpt_clicked(Int_t index_frame);
    void DoAddFile();
    void OnDoubleClick(TGLVEntry*,Int_t);
    void DoFillHistName();
    void hep_clicked();
    void text_clicked();
    Int_t DoFillGraph();
    Int_t DoPlot();
    Int_t DoSave();
    Int_t DoAddEntry();
    Int_t DoRemoveEntry();
    void Multiply_pT(Int_t);

    ClassDef(MyMainFrame, 0)
};
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
class TestFileList {

RQ_OBJECT("TestFileList")

protected:
   TGTransientFrame *fMain;
   TGFileContainer  *fContents;
   TGPopupMenu      *fMenu;
   TGTextButton     *Open;
   TGTextEntry      *Add_file_name;
   const char       *fname;
   TGCompositeFrame *foptions, *ffile, *fpt, *ferrors;

   void DisplayFile(const TString &fname);
   void DisplayDirectory(const TString &fname);
   void DisplayObject(const TString& fname,const TString& name);

public:
   TestFileList(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h);
   virtual ~TestFileList();

   // slots
   void OnDoubleClick(TGLVEntry*,Int_t);
   void DoMenu(Int_t);
   void CloseWindow();
   void DoOpen();
   void hep_clicked();
   void text_clicked();
   void pt_clicked();
   void pt_low_high_clicked();
   void mt_clicked();
   void mt_low_high_clicked();
   void syst_error_clicked();

};
//-----------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------
MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h)
{
    // Constructor
    // Layout
    fL_main = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5 ,1, 1,5);
    fL_framefiles = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 10, 10);
    fL_framebutton = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 10, 10, 10);
    fL_hframe = new TGLayoutHints(kLHintsExpandY, 2,20,10,10);
    fL_Energyframe = new TGLayoutHints(kLHintsExpandX);
    fL_Centralityframe = new TGLayoutHints(kLHintsExpandX,1,1,1,5);
    fL_upperlower = new TGLayoutHints(kLHintsRight,1,1,3);
    fL_button = new TGLayoutHints(kLHintsLeft);

    fMain  = new TGMainFrame(p,w,h); //create main frame

    // Create horizontal composite frame containing Files button
    hframefiles = new TGCompositeFrame(fMain, 60,20, kHorizontalFrame);
    fMain->AddFrame(hframefiles, fL_main);
    Add_file = new TGTextButton(hframefiles, "&Files");
    // When using Files button the local directory opens to browse objects located in file
    Add_file->Connect("Clicked()","MyMainFrame",this,"DoAddFile()");
    hframefiles->AddFrame(Add_file, fL_framefiles);

    // Create composite frame containing text entry and Save button
    hframeoutput = new TGCompositeFrame(fMain, 60,20, kHorizontalFrame);
    fMain->AddFrame(hframeoutput, fL_main);
    // Text entry for output file name
    Add_output_file_name = new TGTextEntry(hframeoutput, "OutputFileName");
    hframeoutput->AddFrame(Add_output_file_name, fL_framefiles);

    Save = new TGTextButton(hframeoutput, "Save");
    // Button Save for saving graphs in output file
    Save->Connect("Clicked()","MyMainFrame",this,"DoSave()");
    hframeoutput->AddFrame(Save, fL_framefiles);

    dNdptpt = new TGCheckButton(hframeoutput, "1/pT*dN/dpT");
    hframeoutput->AddFrame(dNdptpt, fL_framefiles);
    dNdmt = new TGCheckButton(hframeoutput, "dN/dmT");
    hframeoutput->AddFrame(dNdmt, fL_framefiles);
    dNdmtmt = new TGCheckButton(hframeoutput, "1/mT*dN/dmT");
    hframeoutput->AddFrame(dNdmtmt, fL_framefiles);

    ComboPID = new TGComboBox(hframeoutput,200);
    for(Int_t i_pid = 0; i_pid < 22; ++i_pid)
    {
        ComboPID->AddEntry(label_full_pid_spectra[i_pid], i_pid);
    }
    hframeoutput->AddFrame(ComboPID, fL_framefiles);
    ComboPID->Resize(100,20);
    ComboPID->Select(0);

    // Create composite frame containing text buttons: Add entry, Remove entry and Exit
    hframebutton = new TGCompositeFrame(fMain,60,20, kHorizontalFrame);
    fMain->AddFrame(hframebutton,fL_main);

    // Buttons to add or remove entries
    Add_entry = new TGTextButton(hframebutton, "&Add entry");
    Add_entry->Connect("Clicked()","MyMainFrame",this,"DoAddEntry()");
    Remove_entry = new TGTextButton(hframebutton, "&Remove entry");
    Remove_entry->Connect("Clicked()","MyMainFrame",this,"DoRemoveEntry()");
    hframebutton->AddFrame(Add_entry, fL_framebutton);
    hframebutton->AddFrame(Remove_entry, fL_framebutton);

    Plot  = new TGTextButton(hframebutton, "&Plot");
    // When using Plot button the canvas opens
    Plot->Connect("Clicked()","MyMainFrame",this,"DoPlot()");
    hframebutton->AddFrame(Plot, fL_framebutton);

    // Button Exit to close window
    Exit   = new TGTextButton(hframebutton,"&Exit ","gApplication->Terminate()");
    hframebutton->AddFrame(Exit, fL_framebutton);

    fMain->SetWindowName("Select Data");
    fMain->MapSubwindows();
    fMain->Resize(GetDefaultSize());
    fMain->MapWindow();

    h_dummy = new TH1D("h_dummy","h_dummy",400,0,20);    //h_dummy to create a canvas, see function DoPlot
    h_dummy->GetXaxis()->SetRangeUser(-0.03,4.5);
    h_dummy->GetYaxis()->SetRangeUser(-0.1999,0.38);
}
//-------------------------------------------------------------------------------------------------------------------------------------------



//-------------------------------------------------------------------------------------------------------------------------------------------
// Destructor
MyMainFrame::~MyMainFrame()
{

    fMain->Cleanup();
    delete fMain;
}
//-------------------------------------------------------------------------------------------------------------------------------------------



//-------------------------------------------------------------------------------------------------------------------------------------------
void MyMainFrame::DoAddFile()
{
    // Handle Files button
    // Using class TestFileList to browse the local file

    cout<< "DoAddFile"<<endl;
    new TestFileList(gClient->GetRoot(), fMain, 400, 200);
}
//-------------------------------------------------------------------------------------------------------------------------------------------
Int_t MyMainFrame::DoFillGraph()
{
    // Function called if Save or Plot are clicked
    // Fills vector vec_tgae with graphs (TGraphAsymmErrors)

    vec_tgae.clear();
    vec_tgae_syst.clear();
    // text files
    //--------------------------------------------------
    Int_t arr_N_columns[12] = {5,4,4,3, 6, 5, 6,5,5,4,7,6}; // number of columns dependent on the input type

    if(text->IsDown()) // text input file selected
    {
        Int_t i_entry_combo_pid = ComboPID->GetSelected();
        Double_t mass = mass_pid[i_entry_combo_pid];
        cout << "mass: " << mass  << endl;
        TGraphAsymmErrors *tgae = new TGraphAsymmErrors();
        TGraphAsymmErrors *tgae_syst = new TGraphAsymmErrors();

        Int_t input_type = -1;


        if(stat_error->IsDown() && syst_error->IsDown() && !syst_low_high->IsDown() && !pt_error->IsDown() )
        {
            if(pt_high_low ->IsDown() || mt_low_high ->IsDown())  input_type = 0;
            if(pt ->IsDown() || mt ->IsDown())           input_type = 1;
        }
        if(stat_error->IsDown() && syst_error->IsDown() && !syst_low_high->IsDown() && pt_error->IsDown() )
        {
            if(pt_high_low ->IsDown() || mt_low_high ->IsDown())  input_type = 6;
            if(pt ->IsDown() || mt ->IsDown())           input_type = 7;
        }


        if(stat_error->IsDown() && !syst_error->IsDown() && !pt_error->IsDown())
        {
            if(pt_high_low->IsDown() || mt_low_high->IsDown())  input_type = 2;
            if(pt->IsDown() || mt->IsDown())           input_type = 3;
        }
        if(stat_error->IsDown() && !syst_error->IsDown() && pt_error->IsDown())
        {
            if(pt_high_low->IsDown() || mt_low_high->IsDown())  input_type = 8;
            if(pt->IsDown() || mt->IsDown())           input_type = 9;
        }


        if(stat_error->IsDown() && syst_error->IsDown() && syst_low_high->IsDown() && !pt_error->IsDown())
        {
            if(pt_high_low->IsDown() || mt_low_high->IsDown())  input_type = 4;
            if(pt->IsDown() || mt->IsDown())           input_type = 5;
        }
        if(stat_error->IsDown() && syst_error->IsDown() && syst_low_high->IsDown() && pt_error->IsDown())
        {
            if(pt_high_low->IsDown() || mt_low_high->IsDown())  input_type = 10;
            if(pt->IsDown() || mt->IsDown())           input_type = 11;
        }

        if(!(input_type >= 0 && input_type <= 11))
        {
            printf("ERROR: MyMainFrame::DoFillGraph(), wrong input_type!");
            return -1;
        }
        Int_t N_columns = arr_N_columns[input_type];

        //------------------
        // Loop over text file
       
        Int_t i_point = 0;
        Int_t i_block = 0;
        Int_t i_name  = 0;

        while(!myfile.eof())
        {
            std::string str;
            std::getline(myfile,str);

            cout << "Input line: " << str << endl;
            // identifier
            if(str == "# end of file")
            {
                printf("Now at end of file \n");
                vec_tgae.push_back((TGraphAsymmErrors*)tgae->Clone());
                vec_tgae_syst.push_back((TGraphAsymmErrors*)tgae_syst->Clone());
                break;
            }
            std::string sub_str;
            sub_str = str.substr(0,1);
            if(sub_str == "#")
            {
                std::size_t found = str.find("#");
                str.replace(0, found+1, "");

                found = str.find("_");
                sub_str = str.substr(0,found+1);
                HistName = sub_str;
                str.replace(0, found+1, "");

                for (Int_t i_test = 0; i_test < 7; i_test++)
                {
                    found = str.find("_");
                    sub_str = str.substr(0,found+1);
                    HistName += sub_str;
                    str.replace(0, found+1, "");
                }
                sub_str = str.substr(0,1);
                if ( sub_str != " ") HistName+= sub_str;
                str.replace(0, 1, "");
                sub_str = str.substr(0,1);
                if ( sub_str != " ") HistName+= sub_str;
                str.replace(0, 1, "");
                sub_str = str.substr(0,1);
                if ( sub_str != " ") HistName+= sub_str;
                str.replace(0, 1, "");
                sub_str = str.substr(0,1);
                if ( sub_str != " ") HistName+= sub_str;
                //HistName += sub_str;

                //found = str.find("\n");
                //sub_str = str.substr(0,found);
                //HistName = sub_str;
                cout << HistName.Data() << endl;
                if(i_name>=0)
                {
                    vec_HistName.push_back(HistName.Copy());
                    HistName.Clear();
                }
                i_name++;
                i_point = 0;
                if(i_block > 0)
                {
                    printf("Push back new vector A \n");
                    vec_tgae.push_back((TGraphAsymmErrors*)tgae->Clone());
                    vec_tgae_syst.push_back((TGraphAsymmErrors*)tgae_syst->Clone());
                    tgae      ->Set(0);
                    tgae_syst ->Set(0);
                    tgae      ->Clear();
                    tgae_syst ->Clear();
                }
                i_block++;
            }
            else
            {
                vector<Double_t> extracted_values;
                extracted_values = extract_values_from_string(str,N_columns);
                for(Int_t i_value = 0; i_value < (Int_t)extracted_values.size(); i_value++)
                {
                    printf("i_value: %d, value: %4.6f \n",i_value,extracted_values[i_value]);
                }

                Double_t mean           = 0.0;
                Double_t x_low_err      = 0.0;
                Double_t x_up_err       = 0.0;
                Double_t y_low_err_stat = 0.0;
                Double_t y_up_err_stat  = 0.0;
                Double_t y_low_err_syst = 0.0;
                Double_t y_up_err_syst  = 0.0;
                Double_t yvalue         = 0.0;

                if(input_type == 0 || input_type == 2 || input_type == 4   ) // pt high low selected
                {
                    if (pt_high_low->IsDown())
                    {
                        mean = 0.5*(extracted_values[0] + extracted_values[1]);
                        x_low_err = fabs(mean - extracted_values[0]);
                        x_up_err  = fabs(mean - extracted_values[1]);
                    }
                    if (mt_low_high->IsDown())
                    {
                        Double_t x_high = TMath::Sqrt(extracted_values[1]*extracted_values[1]-mass*mass);
                        Double_t x_low = TMath::Sqrt(extracted_values[0]*extracted_values[0]-mass*mass);
                        mean = 0.5*(x_high + x_low);
                        x_low_err = fabs(mean - x_low);
                        x_up_err  = fabs(mean - x_high);
                    }

                    yvalue    = extracted_values[2];
                    y_low_err_stat = extracted_values[3];
                    y_up_err_stat  = extracted_values[3];


                    if(input_type == 0 && !syst_percent->IsDown()) // with systematic error
                    {
                        y_low_err_syst = extracted_values[4];
                        y_up_err_syst  = extracted_values[4];
                    }
                    if(input_type == 0 && syst_percent->IsDown()) // with systematic error percent
                    {
                        y_low_err_syst = extracted_values[4]*yvalue*0.01;
                        y_up_err_syst  = extracted_values[4]*yvalue*0.01;
                    }
                    if(input_type == 4) // with systematic error low_high
                    {
                        y_low_err_syst = extracted_values[4];
                        y_up_err_syst  = extracted_values[5];
                    }
                    if(input_type == 2 || input_type == 8) // without systematic error
                    {
                        y_low_err_syst = 0.0;
                        y_up_err_syst  = 0.0;
                    }
                }
                else // pt value direct
                {
                    if (pt->IsDown()) mean    = extracted_values[0];
                    if (mt->IsDown()) mean    = TMath::Sqrt(extracted_values[0]*extracted_values[0]-mass*mass);
                    if (!pt_error->IsDown())
                    {
                        yvalue  = extracted_values[1];
                        y_low_err_stat = extracted_values[2];
                        y_up_err_stat  = extracted_values[2];
                    }

                    if (pt_error->IsDown())
                    {
                        printf("Test A \n");
                        yvalue  = extracted_values[2];
                        y_low_err_stat = extracted_values[3];
                        y_up_err_stat  = extracted_values[3];
                        x_low_err = extracted_values[1];
                        x_up_err  = extracted_values[1];

                    }
                    if (mt->IsDown() && pt_error->IsDown())
                    {
                        x_low_err = TMath::Sqrt(extracted_values[1]*extracted_values[1]-mass*mass);
                        x_up_err  = TMath::Sqrt(extracted_values[1]*extracted_values[1]-mass*mass);


                    }

                    if(input_type == 1 && !syst_percent->IsDown()) // with systematic error
                    {
                        printf("Test B \n");
                        y_low_err_syst = extracted_values[3];
                        y_up_err_syst  = extracted_values[3];

                    }
                    if(input_type == 7 && !syst_percent->IsDown()) // with systematic error
                    {
                        printf("Test C \n");
                        y_low_err_syst = extracted_values[4];
                        y_up_err_syst  = extracted_values[4];
                    }
                    if(input_type == 5 && !syst_percent->IsDown()) // with systematic error
                    {
                        printf("Test D \n");
                        y_low_err_syst = extracted_values[3];
                        y_up_err_syst  = extracted_values[4];
                    }
                    if(input_type == 11 && !syst_percent->IsDown()) // with systematic error
                    {
                        printf("Test E \n");
                        y_low_err_syst = extracted_values[4];
                        y_up_err_syst  = extracted_values[5];
                    }
                    if(input_type == 1 && syst_percent->IsDown()) // with systematic error percent
                    {
                        printf("Test F \n");
                        y_low_err_syst = extracted_values[3]*yvalue*0.01;
                        y_up_err_syst  = extracted_values[3]*yvalue*0.01;

                    }
                    if(input_type == 7 && syst_percent->IsDown()) // with systematic error percent
                    {
                        printf("Test G \n");
                        y_low_err_syst = extracted_values[4]*yvalue*0.01;
                        y_up_err_syst  = extracted_values[4]*yvalue*0.01;
                    }
                    if(input_type == 3 || input_type == 9) // without systematic error
                    {
                        printf("Test H \n");
                        y_low_err_syst = 0.0;
                        y_up_err_syst  = 0.0;
                    }

                }

                if (dNdmt->IsDown())
                {
                    printf("Test A \n");
                    yvalue         *= mean*1/(TMath::Sqrt(mean*mean+mass*mass)); // dNdmT
                    y_low_err_stat *= mean*1/(TMath::Sqrt(mean*mean+mass*mass));
                    y_up_err_stat  *= mean*1/(TMath::Sqrt(mean*mean+mass*mass));
                    y_low_err_syst *= mean*1/(TMath::Sqrt(mean*mean+mass*mass));
                    y_up_err_syst  *= mean*1/(TMath::Sqrt(mean*mean+mass*mass));
                    //cout << "yvalue " << yvalue<< "y_low_error_stat " << y_low_err_stat<< "y_up_err_stat " << y_up_err_stat << endl;
                }
                if (dNdmtmt->IsDown())
                {
                    printf("Test B \n");
                    yvalue *= TMath::Sqrt(mean*mean+mass*mass); // * mT
                    y_low_err_stat *=  TMath::Sqrt(mean*mean+mass*mass) ;
                    y_up_err_stat  *=  TMath::Sqrt(mean*mean+mass*mass);
                    y_low_err_syst *=  TMath::Sqrt(mean*mean+mass*mass);
                    y_up_err_syst  *=  TMath::Sqrt(mean*mean+mass*mass);
                    //cout << "yvalue " << yvalue<< "y_low_error_stat " << y_low_err_stat<< "y_up_err_stat " << y_up_err_stat << endl;
                }
                if (dNdptpt->IsDown())
                {
                    yvalue *= mean; // * pT
                    y_low_err_stat *= mean;
                    y_up_err_stat *= mean;
                    y_low_err_syst *= mean;
                    y_up_err_syst *= mean;
                }

                tgae      ->SetPoint(i_point, mean, yvalue);
                tgae_syst ->SetPoint(i_point, mean, yvalue);
                tgae      ->SetPointError(i_point, x_low_err, x_up_err, y_low_err_stat, y_up_err_stat);
                tgae_syst ->SetPointError(i_point, x_low_err, x_up_err, y_low_err_syst, y_up_err_syst);
                i_point++;

                cout << "xvalue: " << mean << ", x_low_error: " << x_low_err<< ", x_up_err: " << x_up_err << ", yvalue: " << yvalue<< ", y_low_error_stat: " << y_low_err_stat<< ", y_up_err_stat: " << y_up_err_stat << endl;

            }

        }
        //------------------

    }
    // End of text input file selected
    //--------------------------------------------------

#if 0
    if(text->IsDown()) // text input file selected
    {
        Int_t N_columns = 0;

        if(stat_error->IsDown()&& !syst_error->IsDown())
        {
            TGraphAsymmErrors *tgae = new TGraphAsymmErrors();

            Int_t i_point = 0;
            Int_t i_block = 0;
            Int_t i_name = 0;

            while(!myfile.eof())
            {
                std::string str;
                std::getline(myfile,str);
                cout << "Line A: " << str << endl;
                // identifier
                if(str == "# end of file")
                {
                    printf("Now at end of file \n");
                    break;
                }
                std::string sub_str;
                sub_str = str.substr(0,1);
                if(sub_str == "#")
                {
                    std::size_t found = str.find("#");
                    str.replace(0, found+1, "");
                    found = str.find("\n");
                    sub_str = str.substr(0,found);
                    HistName = sub_str;
                    cout << HistName.Data() << endl;
                    if(i_name>=0)
                    {
                        vec_HistName.push_back(HistName.Copy());
                        HistName.Clear();
                    }
                    i_name++;
                    i_point = 0;
                    if(i_block > 0)
                    {
                        printf("Push back new vector A \n");
                        vec_tgae.push_back((TGraphAsymmErrors*)tgae->Clone());
                        tgae ->Set(0);
                        tgae ->Clear();
                    }
                    i_block++;
                }
                else
                {
                    if(pt_high_low->IsDown())
                    {
                        Double_t values[4];

                        std::size_t found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        sub_str = str.substr(0,found);
                        values[0] = std::stod(sub_str);

                        str.replace(0,found+1,"");
                        foundB = str.find_first_not_of("\t");
                        foundC = str.find_first_not_of(" ");
                        if(foundC > foundB) foundB = foundC;
                        if(foundB > 0) str.replace(0,foundB,"");

                        found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        sub_str = str.substr(0,found);
                        values[1] = std::stod(sub_str);

                        str.replace(0,found+1,"");
                        foundB = str.find_first_not_of("\t");
                        foundC = str.find_first_not_of(" ");
                        if(foundC > foundB) foundB = foundC;
                        if(foundB > 0) str.replace(0,foundB,"");

                        found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        sub_str = str.substr(0,found);
                        values[2] = std::stod(sub_str);

                        str.replace(0,found+1,"");
                        foundB = str.find_first_not_of("\t");
                        foundC = str.find_first_not_of(" ");
                        if(foundC > foundB) foundB = foundC;
                        if(foundB > 0) str.replace(0,foundB,"");

                        found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        std::size_t found_end = str.length();
                        sub_str = str.substr(0,found_end);
                        values[3] = std::stod(sub_str);

                        Double_t yvalue;
                        Double_t mean = 0.5*(values[0]+values[1]);
                        Double_t low_err = fabs(mean - values[0]);
                        Double_t up_err = fabs(mean - values[1]);

                        if ( dNdptpt->IsDown())
                        {
                            yvalue = values[2]*mean;
                        }
                        else
                        {
                            yvalue = values[2];

                        }
                        tgae->SetPoint(i_point, mean, yvalue);
                        tgae->SetPointError(i_point, low_err,up_err,values[3],values[3]);
                        i_point++;
                        printf("values: {%4.5f, %4.5f, %4.5f, %4.5f} \n",values[0],values[1],values[2],values[3]);
                    }
                    if (pt->IsDown())
                    {
                        //printf("pt is down \n");
                        Double_t values[3];

                        std::size_t found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        sub_str = str.substr(0,found);
                        cout << "A1 sub_str: " << sub_str << ", found: " <<  found << ", foundC: " << foundC << endl;
                        values[0] = std::stod(sub_str);

                        str.replace(0,found+1,"");
                        foundB = str.find_first_not_of("\t");
                        foundC = str.find_first_not_of(" ");
                        cout << str << ", foundB: " << foundB << ", foundC: " << foundC << endl;
                        if(foundC > foundB) foundB = foundC;
                        if(foundB > 0) str.replace(0,foundB,"");

                        found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        sub_str = str.substr(0,found);
                        cout << "B1 sub_str: " << sub_str << ", found: " <<  found << ", foundC: " << foundC << endl;
                        values[1] = std::stod(sub_str);

                        str.replace(0,found+1,"");
                        foundB = str.find_first_not_of("\t");
                        foundC = str.find_first_not_of(" ");
                        if(foundC > foundB) foundB = foundC;
                        if(foundB > 0) str.replace(0,foundB,"");

                        found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        std::size_t found_end = str.length();
                        sub_str = str.substr(0,found_end);
                        cout << "C1 sub_str: " << sub_str << ", found: " <<  found << ", foundC: " << foundC << endl;
                        values[2] = std::stod(sub_str);

                        Double_t yvalue;
                        if ( dNdptpt->IsDown())
                        {
                            yvalue = values[1]*values[0];
                        }
                        else
                        {
                            yvalue = values[1];

                        }
                        tgae->SetPoint(i_point, values[0], yvalue);
                        tgae->SetPointError(i_point, 0.0, 0.0, values[2], values[2]);
                        i_point++;
                        printf("values: {%4.5f, %4.5f, %4.5f} \n",values[0],values[1],values[2]);
                    }

                }
            }
        }




        if (stat_error->IsDown() && syst_error->IsDown())
        {

            TGraphAsymmErrors *tgae = new TGraphAsymmErrors();
            TGraphAsymmErrors *tgae_syst = new TGraphAsymmErrors();

            Int_t i_point = 0;
            Int_t i_block = 0;
            Int_t i_name = 0;

            while(!myfile.eof())
            {
                std::string str;
                std::getline(myfile,str);
                // identifier
                if(str == "# end of file")
                {
                    printf("Now at end of file \n");
                    break;
                }
                std::string sub_str;
                sub_str = str.substr(0,1);
                if(sub_str == "#")
                {
                    std::size_t found = str.find("#");
                    str.replace(0, found+1, "");
                    found = str.find("\n");
                    sub_str = str.substr(0,found);
                    HistName = sub_str;
                    cout << HistName.Data() << endl;
                    if(i_name>=0)
                    {
                        vec_HistName.push_back(HistName.Copy());
                        HistName.Clear();
                    }
                    i_name++;
                    i_point = 0;
                    if(i_block > 0)
                    {
                        printf("Push back new vector \n");
                        vec_tgae.push_back((TGraphAsymmErrors*)tgae->Clone());
                        vec_tgae_syst.push_back((TGraphAsymmErrors*)tgae_syst->Clone());
                        tgae ->Set(0);
                        tgae_syst ->Set(0);
                        tgae ->Clear();
                        tgae_syst ->Clear();
                    }
                    i_block++;
                }

                else
                {

                    if (pt_high_low->IsDown())
                    {
                        cout << "A: " << str << endl;

                        Double_t values[5];
                        std::size_t found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        sub_str = str.substr(0,found);
                        cout << "A1 sub_str: " << sub_str << ", found: " <<  found << ", foundC: " << foundC << endl;
                        values[0] = std::stod(sub_str);

                        str.replace(0,found+1,"");
                        cout << "A2: " << str << endl;
                        foundB = str.find_first_not_of("\t");
                        foundC = str.find_first_not_of(" ");
                        if(foundC > foundB) foundB = foundC;
                        if(foundB > 0) str.replace(0,foundB,"");

                        cout << "B: " << str << endl;

                        found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        sub_str = str.substr(0,found);
                        cout << "B1 sub_str: " << sub_str << ", found: " <<  found << ", foundC: " << foundC << endl;
                        values[1] = std::stod(sub_str);

                        str.replace(0,found+1,"");
                        foundB = str.find_first_not_of("\t");
                        foundC = str.find_first_not_of(" ");
                        if(foundC > foundB) foundB = foundC;
                        if(foundB > 0) str.replace(0,foundB,"");

                        cout << "C: " << str << endl;

                        found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        sub_str = str.substr(0,found);
                        cout << "C1 sub_str: " << sub_str << ", found: " <<  found << ", foundC: " << foundC << endl;
                        values[2] = std::stod(sub_str);

                        str.replace(0,found+1,"");
                        foundB = str.find_first_not_of("\t");
                        foundC = str.find_first_not_of(" ");
                        if(foundC > foundB) foundB = foundC;
                        if(foundB > 0) str.replace(0,foundB,"");
                        found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        sub_str = str.substr(0,found);
                        values[3] = std::stod(sub_str);

                        str.replace(0,found+1,"");
                        foundB = str.find_first_not_of("\t");
                        foundC = str.find_first_not_of(" ");
                        if(foundC > foundB) foundB = foundC;
                        if(foundB > 0) str.replace(0,foundB,"");

                        found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        std::size_t found_end = str.length();
                        sub_str = str.substr(0,found_end);
                        values[4] = std::stod(sub_str);

                        Double_t yvalue;
                        Double_t mean = 0.5*(values[0]+values[1]);
                        Double_t low_err = fabs(mean - values[0]);
                        Double_t up_err = fabs(mean - values[1]);

                        if ( dNdptpt->IsDown())
                        {
                            yvalue = values[2]*mean;
                        }
                        else
                        {
                            yvalue = values[2];

                        }
                        tgae->SetPoint(i_point, mean, yvalue);
                        tgae_syst->SetPoint(i_point, mean, yvalue);
                        tgae->SetPointError(i_point, low_err,up_err,values[3],values[3]);
                        tgae_syst->SetPointError(i_point, low_err,up_err, values[4], values[4]);
                        //tgae->SetPointError(i_point, 0.0,0.0,values[3],values[3]);
                        //tgae_syst->SetPointError(i_point, 0.0,0.0, values[4], values[4]);
                        i_point++;
                        printf("values: {%4.5f, %4.5f, %4.5f, %4.5f. %4.5f} \n",values[0],values[1],values[2],values[3], values[4]);
                    }
                    if (pt->IsDown())
                    {
                        Double_t values[4];

                        std::size_t found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        sub_str = str.substr(0,found);
                        values[0] = std::stod(sub_str);

                        str.replace(0,found+1,"");
                        foundB = str.find_first_not_of("\t");
                        foundC = str.find_first_not_of(" ");
                        if(foundC > foundB) foundB = foundC;
                        if(foundB > 0) str.replace(0,foundB,"");

                        found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        sub_str = str.substr(0,found);
                        values[1] = std::stod(sub_str);

                        str.replace(0,found+1,"");
                        foundB = str.find_first_not_of("\t");
                        foundC = str.find_first_not_of(" ");
                        if(foundC > foundB) foundB = foundC;
                        if(foundB > 0) str.replace(0,foundB,"");

                        found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        sub_str = str.substr(0,found);
                        values[2] = std::stod(sub_str);

                        str.replace(0,found+1,"");
                        foundB = str.find_first_not_of("\t");
                        foundC = str.find_first_not_of(" ");
                        if(foundC > foundB) foundB = foundC;
                        if(foundB > 0) str.replace(0,foundB,"");

                        found = str.find("\t");
                        foundC = str.find(" ");
                        if(foundC < found) found = foundC;
                        std::size_t found_end = str.length();
                        sub_str = str.substr(0,found_end);
                        values[3] = std::stod(sub_str);

                        Double_t yvalue;
                        if ( dNdptpt->IsDown())
                        {
                            yvalue = values[1]*values[0];
                        }
                        else
                        {
                            yvalue = values[1];

                        }
                        tgae->SetPoint(i_point, values[0], yvalue);
                        tgae_syst->SetPoint(i_point, values[0], yvalue);
                        tgae->SetPointError(i_point, 0.0, 0.0, values[2], values[2]);
                        tgae_syst->SetPointError(i_point, 0.0, 0.0, values[3], values[3]);
                        i_point++;
                        printf("values: {%4.5f, %4.5f, %4.5f, %4.5f} \n",values[0],values[1],values[2],values[3]);
                    }
                }
            }
        }
    }
#endif
   // printf("Test E \n");
    if (hep->IsDown())
    {
     //   printf("Test F \n");
        for(Int_t i_entry=0 ; i_entry < (Int_t)vec_fCombo.size(); i_entry++)
        {
            TString table_name;
            Int_t i_select_table = vec_table[i_entry]->GetNumber();   //Get number of selected table
            if(i_select_table == 0) return 0;
            table_name += "Table ";
            table_name += table_name.Format("%d", i_select_table);
            table_name += "/Graph1D_y";
            Int_t i_select_graph = vec_graph[i_entry]->GetNumber();   //Get number of selected graph
            table_name += table_name.Format("%d", i_select_graph);
            cout <<table_name <<endl;
            vec_tgae.push_back( (TGraphAsymmErrors*)newfile->Get(table_name.Data())); //Fill vector vec_tgae
            if(dNdptpt->IsDown())
            {
                Multiply_pT(i_entry);

            }
           // for (Int_t i_point = 0; i_point < vec_tgae[i_entry] ->GetN(); i_point++)
           // {

             //   vec_tgae[i_entry]->GetPoint(i_point,pT,y_val);
        //}
        }
    }
    return 1;

}
//-------------------------------------------------------------------------------------------------------------------------------------------
void MyMainFrame::Multiply_pT(Int_t i_entry)
{
    for (Int_t i_point = 0; i_point < vec_tgae[i_entry] ->GetN(); i_point++)
    {
        //printf("Test E \n");
        Double_t pT,y_val,err_X_high,err_X_low,err_Y_high,err_Y_low;
        vec_tgae[i_entry]->GetPoint(i_point,pT,y_val);
        err_X_high = fabs(vec_tgae[i_entry] ->GetErrorXhigh(i_point));
        err_X_low  = fabs(vec_tgae[i_entry] ->GetErrorXlow(i_point));
        err_Y_high = fabs(vec_tgae[i_entry] ->GetErrorYhigh(i_point));
        err_Y_low  = fabs(vec_tgae[i_entry] ->GetErrorYlow(i_point));
        vec_tgae[i_entry] ->SetPoint(i_point,pT,y_val*pT);
        vec_tgae[i_entry] ->SetPointEYhigh(i_point,err_Y_high*pT);
        vec_tgae[i_entry] ->SetPointEYlow(i_point,err_Y_low*pT);
        vec_tgae[i_entry] ->SetPointEXhigh(i_point,err_X_high);
        vec_tgae[i_entry] ->SetPointEXlow(i_point,err_X_low);
    }


}
//-------------------------------------------------------------------------------------------------------------------------------------------
void MyMainFrame::DoFillHistName()
{
    // Function called if Save is clicked
    // Creates names for every graph

    vec_HistName.clear();
    for(Int_t i_entry=0 ; i_entry < (Int_t)vec_fCombo.size(); i_entry++)
    {
        if (vec_v2[i_entry]->IsDown()) HistName = "v2_";
        if (vec_dNdpt[i_entry]->IsDown()) HistName = "dNdpt_";

        Int_t i_select = vec_fCombo[i_entry]->GetSelected();
        HistName += "PID_";
        if (i_select == 0)  HistName += "Pi+";
        if (i_select == 1)  HistName += "Pi-";
        if (i_select == 2)  HistName += "K+";
        if (i_select == 3)  HistName += "K-";
        if (i_select == 4)  HistName += "P";
        if (i_select == 5)  HistName += "Pbar";
        if (i_select == 6)  HistName += "Phi";
        if (i_select == 7)  HistName += "Lambda";
        if (i_select == 8)  HistName += "Lambdabar";
        if (i_select == 9)  HistName += "Xi-";
        if (i_select == 10) HistName += "Xibar+";
        if (i_select == 11) HistName += "Omega-";
        if (i_select == 12) HistName += "Omegabar+";
        if (i_select == 13) HistName += "D0";
        if (i_select == 14) HistName += "J/Psi";
        if (i_select == 15) HistName += "K0S";
        if (i_select == 16) HistName += "d";
        if (i_select == 17) HistName += "dbar";
        if (i_select == 18) HistName += "He3";
        if (i_select == 19) HistName += "He3bar";
        if (i_select == 20) HistName += "t";
        if (i_select == 21) HistName += "Upsilon";
        if (i_select == 22) HistName += "LambdaC";

        Double_t i_select_energy = vec_energy_entry[i_entry]->GetNumber();
        HistName += "_E_";
        HistName += HistName.Format("%.1f", i_select_energy);

        Double_t i_select_C_lower = vec_centrality_lower[i_entry]->GetNumber();
        HistName += "_C_";
        HistName += i_select_C_lower;
        //HistName += HistName.Format("%d", i_select_C_lower);

        Double_t i_select_C_upper = vec_centrality_upper[i_entry]->GetNumber();
        HistName += "_";
        HistName += i_select_C_upper;
         HistName += "_stat";
        //HistName += HistName.Format("%d", i_select_C_upper);

        vec_HistName.push_back(HistName.Copy());           //Fill vector vec_HistName
        cout << HistName<<endl;
        HistName.Clear();
    }
}
//-------------------------------------------------------------------------------------------------------------------------------------------
Int_t MyMainFrame::DoSave()
{
    // Handle Save button
    // Saves the graphs with selected names in root file

    TString output_file_name = Add_output_file_name->GetDisplayText();          //Name of output file = name displayed in text entry
    output_file_name += ".root";
    cout << "output_file_name: " << output_file_name.Data() << endl;
    RootFile = new TFile(output_file_name.Data(), "RECREATE");          //Create new root file
    if ( RootFile->IsOpen() ) printf("File opened successfully\n");
    DoFillGraph();
    cout << vec_tgae.size() << endl;
    cout << vec_HistName.size() << endl;

    if (hep->IsDown())
    {
        DoFillHistName();
        if (N_used_entries == 0) return 0;    // Add entry first
        for(Int_t i_entry=0 ; i_entry < N_used_entries; i_entry++)
        {
            Int_t i_select_table = vec_table[i_entry]->GetNumber();
            if(i_select_table == 0) return 0;                               // Select table before saving

            RootFile->cd();
            vec_tgae[i_entry]->Write(vec_HistName[i_entry].Data());         // Save graphs in root file
        }
        printf("DoSave, number of entries: %d \n",N_used_entries);
    }


    if (text->IsDown())
    {
        //if (stat_error->IsDown() && syst_error->IsDown())
        //{
            for(Int_t index_vector = 0; index_vector < (Int_t)vec_tgae.size(); index_vector++)
            {
                RootFile->cd();
                TString HistName_stat = vec_HistName[index_vector];
                HistName_stat += "_stat";
                TString HistName_syst = vec_HistName[index_vector];
                HistName_syst += "_syst";
                vec_tgae[index_vector]->Write(HistName_stat);
                vec_tgae_syst[index_vector]->Write(HistName_syst);
                cout  << HistName_stat.Data() << endl;
                cout  << HistName_syst.Data() << endl;
            }

       // }
       /*
        if (stat_error->IsDown() && !syst_error->IsDown())
        {
            for(Int_t index_vector = 0; index_vector < (Int_t)vec_tgae.size(); index_vector++)
            {
                RootFile->cd();
                TString HistName_stat = vec_HistName[index_vector];
                HistName_stat += "_stat";
                vec_tgae[index_vector]->Write(HistName_stat);
            }

        } */
    }

    return 1;
}
//-------------------------------------------------------------------------------------------------------------------------------------------
void MyMainFrame::CloseWindow()
{
    delete this;
}
//-------------------------------------------------------------------------------------------------------------------------------------------
void MyMainFrame::v2_clicked(Int_t index_frame)
{
    // Handle v2 button
    // Needed because only v2 or dNdpt should be marked

    printf("v2_clicked \n");
    vec_dNdpt[index_frame] ->SetState(kButtonUp);
}
//-------------------------------------------------------------------------------------------------------------------------------------------
void MyMainFrame::dNdpt_clicked(Int_t index_frame)
{
    // Handle dNdpt button
    // See v2_clicked()

    printf("dNdpt_clicked \n");
    vec_v2[index_frame] ->SetState(kButtonUp);
}
//-------------------------------------------------------------------------------------------------------------------------------------------
void MyMainFrame::hep_clicked()
{
    // Handle hep file button
    // See v2_clicked()

    printf("hep_clicked \n");
    text ->SetState(kButtonUp);
}
//-------------------------------------------------------------------------------------------------------------------------------------------
void MyMainFrame::text_clicked()
{
    // Handle text file button
    // See v2_clicked()

    printf("text_clicked \n");
    hep->SetState(kButtonUp);
}
//-------------------------------------------------------------------------------------------------------------------------------------------
Int_t MyMainFrame::DoAddEntry()
{
    // Handle Add entry button
    // Creates buttons and number entries to select names and table numbers
    if (text->IsDown()) return 0;
    if (N_used_entries > 8) return 0;        //Max. 9 entries are possible

    if(!flag_file_opened) return 0;

    printf("Do_Add_entry \n");
    Int_t index_frame =  (Int_t)vec_hframe.size();
    if(index_frame <= N_used_entries)
    {
        vec_hframe.push_back(new TGCompositeFrame(fMain,60,20, kHorizontalFrame));
        vec_V2frame.push_back(new TGCompositeFrame(vec_hframe[index_frame], kVerticalFrame));
        vec_PIDframe.push_back(new TGCompositeFrame(vec_hframe[index_frame], kVerticalFrame));
        vec_Tableframe.push_back(new TGCompositeFrame(vec_hframe[index_frame], kVerticalFrame));
        vec_Energyframe.push_back(new TGCompositeFrame(vec_hframe[index_frame], kVerticalFrame));
        vec_Centrality_labelframe.push_back(new TGCompositeFrame(vec_hframe[index_frame], kVerticalFrame));
        vec_Centralityframe.push_back(new TGCompositeFrame(vec_hframe[index_frame], kVerticalFrame));
        vec_Graphframe.push_back(new TGCompositeFrame(vec_hframe[index_frame], kVerticalFrame));

        // Buttons to select v2 or dNdpt
        vec_fButtonGroup.push_back(new TGVButtonGroup(vec_V2frame[index_frame]));
        vec_v2.push_back(new TGRadioButton(vec_V2frame[index_frame], "v2"));
        vec_dNdpt.push_back( new TGRadioButton(vec_V2frame[index_frame], "dN/dpt"));

        // Number entry (integer) to select table and graph number and centrality
        vec_table.push_back(new TGNumberEntry(vec_Tableframe[index_frame],0,9,99, TGNumberFormat::kNESInteger,
                                              TGNumberFormat::kNEANonNegative,
                                              TGNumberFormat::kNELLimitMinMax,
                                              0, 99999));
        vec_graph.push_back(new TGNumberEntry(vec_Graphframe[index_frame],1,9,99, TGNumberFormat::kNESInteger,
                                              TGNumberFormat::kNEANonNegative,
                                              TGNumberFormat::kNELLimitMinMax,
                                              0, 99999));
        vec_centrality_lower.push_back(new TGNumberEntry(vec_Centralityframe[index_frame],0,9,99, TGNumberFormat::kNESInteger,
                                              TGNumberFormat::kNEANonNegative,
                                              TGNumberFormat::kNELLimitMinMax,
                                              0, 99999));
        vec_centrality_upper.push_back(new TGNumberEntry(vec_Centralityframe[index_frame],0,9,99, TGNumberFormat::kNESInteger,
                                              TGNumberFormat::kNEANonNegative,
                                              TGNumberFormat::kNELLimitMinMax,
                                              0, 99999));

        // Number(float) entries to select energy
        vec_energy_entry.push_back(new TGNumberEntry(vec_Energyframe[index_frame],0.1,TGNumberFormat::kNESRealOne));

        // Combobox for PID
        vec_fCombo.push_back(new TGComboBox(vec_PIDframe[index_frame]));

        // Labels
        vec_label_PID.push_back(new TGLabel(vec_PIDframe[index_frame], "PID"));
        vec_label_Energy.push_back(new TGLabel(vec_Energyframe[index_frame], "Energy"));
        vec_label_Table.push_back(new TGLabel(vec_Tableframe[index_frame], "Table"));
        vec_label_Graph.push_back(new TGLabel(vec_Graphframe[index_frame], "Graph"));
        vec_label_Centrality.push_back(new TGLabel(vec_Centralityframe[index_frame], "Centrality"));
        vec_label_upper.push_back(new TGLabel(vec_Centrality_labelframe[index_frame], "upper"));
        vec_label_lower.push_back(new TGLabel(vec_Centrality_labelframe[index_frame], "lower"));
        vec_label_space.push_back(new TGLabel(vec_Centrality_labelframe[index_frame], ""));


        // Add entries and buttons to frames
        vec_fButtonGroup[index_frame]->AddFrame( vec_v2[index_frame]);
        vec_fButtonGroup[index_frame]->AddFrame( vec_dNdpt[index_frame]);
        vec_v2[index_frame]->Connect("Clicked()", "MyMainFrame", this, Form("v2_clicked(Int_t = %d)", index_frame));
        vec_dNdpt[index_frame]->Connect("Clicked()", "MyMainFrame", this, Form("dNdpt_clicked(Int_t = %d)", index_frame));
        vec_V2frame[index_frame]->AddFrame(vec_fButtonGroup[index_frame]);
        vec_PIDframe[index_frame]->AddFrame(vec_label_PID[index_frame]);
        vec_PIDframe[index_frame]->AddFrame(vec_fCombo[index_frame]);
        vec_fCombo[index_frame]->AddEntry("Pi+", 0);
        vec_fCombo[index_frame]->AddEntry("Pi-", 1);
        vec_fCombo[index_frame]->AddEntry("K+", 2);
        vec_fCombo[index_frame]->AddEntry("K-", 3);
        vec_fCombo[index_frame]->AddEntry("P", 4);
        vec_fCombo[index_frame]->AddEntry("Pbar", 5);
        vec_fCombo[index_frame]->AddEntry("Phi", 6);
        vec_fCombo[index_frame]->AddEntry("Lambda", 7);
        vec_fCombo[index_frame]->AddEntry("Lambdabar", 8);
        vec_fCombo[index_frame]->AddEntry("Xi-", 9);
        vec_fCombo[index_frame]->AddEntry("Xibar+", 10);
        vec_fCombo[index_frame]->AddEntry("Omega-", 11);
        vec_fCombo[index_frame]->AddEntry("Omegabar+", 12);
        vec_fCombo[index_frame]->AddEntry("D0", 13);
        vec_fCombo[index_frame]->AddEntry("J/Psi", 14);
        vec_fCombo[index_frame]->AddEntry("K0S", 15);
        vec_fCombo[index_frame]->AddEntry("d", 16);
        vec_fCombo[index_frame]->AddEntry("dbar", 17);
        vec_fCombo[index_frame]->AddEntry("He3", 18);
        vec_fCombo[index_frame]->AddEntry("He3bar", 19);
        vec_fCombo[index_frame]->AddEntry("t", 20);
        vec_fCombo[index_frame]->AddEntry("Upsilon", 21);
        vec_fCombo[index_frame]->AddEntry("LambdaC", 22);
        vec_fCombo[index_frame]->Resize(100, 20);

        vec_Tableframe[index_frame]->AddFrame(vec_label_Table[index_frame]);
        vec_Tableframe[index_frame]->AddFrame(vec_table[index_frame], fL_Energyframe);

        vec_Graphframe[index_frame]->AddFrame(vec_label_Graph[index_frame]);
        vec_Graphframe[index_frame]->AddFrame(vec_graph[index_frame], fL_Energyframe);

        vec_Energyframe[index_frame]->AddFrame(vec_label_Energy[index_frame]);
        vec_Energyframe[index_frame]->AddFrame(vec_energy_entry[index_frame], fL_Energyframe);
        vec_Energyframe[index_frame]->Resize(100,50);

        vec_Centrality_labelframe[index_frame]->AddFrame(vec_label_upper[index_frame], fL_upperlower);
        vec_Centrality_labelframe[index_frame]->AddFrame(vec_label_lower[index_frame], fL_upperlower);

        vec_Centralityframe[index_frame]->AddFrame(vec_label_Centrality[index_frame]);
        vec_Centralityframe[index_frame]->AddFrame(vec_centrality_upper[index_frame], fL_Centralityframe);
        vec_Centralityframe[index_frame]->AddFrame(vec_centrality_lower[index_frame],fL_Centralityframe );

        vec_hframe[index_frame]->AddFrame(vec_V2frame[index_frame],fL_hframe );
        vec_hframe[index_frame]->AddFrame(vec_PIDframe[index_frame],fL_hframe );
        vec_hframe[index_frame]->AddFrame(vec_Tableframe[index_frame],fL_hframe );
        vec_hframe[index_frame]->AddFrame(vec_Graphframe[index_frame],fL_hframe );
        vec_hframe[index_frame]->AddFrame(vec_Energyframe[index_frame],fL_hframe);
        vec_hframe[index_frame]->AddFrame(vec_Centrality_labelframe[index_frame], new TGLayoutHints( 0, 0, 0, 31));
        vec_hframe[index_frame]->AddFrame(vec_Centralityframe[index_frame],fL_hframe);

        vec_Centralityframe[index_frame]->Resize(100, 50);
        vec_Tableframe[index_frame]->Resize(100,50);
        vec_Graphframe[index_frame]->Resize(100,50);
        vec_hframe[index_frame]->Resize(650,200);

        fMain->AddFrame(vec_hframe[index_frame]);
        fMain->Resize();
        fMain->MapSubwindows();
        fMain->MapWindow();
    }

    // need no new buttons and number enrties if they have already been created and removed
    else
    {
        fMain->AddFrame(vec_hframe[N_used_entries]);
        fMain->Resize();
        fMain->MapSubwindows();
        fMain->MapWindow();
    }
    N_used_entries++;
    return 1;
}
//-------------------------------------------------------------------------------------------------------------------------------------------
Int_t MyMainFrame::DoRemoveEntry()
{
    // Handle Remove entry button

    if (N_used_entries==0) return 0;
    printf("Do_Remove_entry \n");
    fMain->RemoveFrame(vec_hframe[N_used_entries-1]);
    fMain->Resize(GetDefaultSize());
    fMain->MapSubwindows();
    fMain->MapWindow();
    N_used_entries--;


    return 1;
}
//-------------------------------------------------------------------------------------------------------------------------------------------
Int_t MyMainFrame::DoPlot()
{
    // Handle Plot button
    // Create new canvas 3X3 and plot the graphs

    if(N_used_entries ==0) return 0;    // Add entry first

    // create new canvas only first time the button is clicked
    if(!c_3X3)
    {
        c_3X3 = new TCanvas("c_3X3","c_3X3",100,200,1400,900);
        c_3X3->SetTopMargin(0.15);
        c_3X3->SetBottomMargin(0.25);
        c_3X3->SetRightMargin(0.15);
        c_3X3->SetLeftMargin(0.25);
        c_3X3->SetLogy(0);
        c_3X3->Divide(3,3,0.0,0.0); // x divide, y divide, x margin, y margin
    }

    // Scaling different pads and plot h_dummy
    for(Int_t i_pad = 0; i_pad <  9; i_pad++)
    {
        Double_t scaling_factor_3X3 = 1.0;
        Double_t offset_3X3 = 0.015;
        Double_t Label_size_3X3     = 0.1;
        if(i_pad == 6)
        {
            scaling_factor_3X3 = 0.78;
            offset_3X3 = 0.015;
        }

        c_3X3->cd(i_pad+1)->SetTicks(1,1);
        c_3X3->cd(i_pad+1)->SetGrid(0,0);
        c_3X3->cd(i_pad+1)->SetFillColor(10);
        c_3X3->cd(i_pad+1)->SetRightMargin(0.01);
        c_3X3->cd(i_pad+1)->SetTopMargin(0.01);
        if(i_pad < 6) c_3X3->cd(i_pad -8)->SetBottomMargin(0.0045);
        h_dummy->SetTitle("");
        h_dummy->SetStats(0);
        h_dummy->GetXaxis()->SetTitleOffset(0.9/scaling_factor_3X3);
        h_dummy->GetYaxis()->SetTitleOffset(1/scaling_factor_3X3);
        h_dummy->GetXaxis()->SetLabelOffset(0.0*scaling_factor_3X3 + offset_3X3);
        if(i_pad == 0) h_dummy->GetXaxis()->SetLabelOffset(0.0*scaling_factor_3X3);
        h_dummy->GetYaxis()->SetLabelOffset(0.01*scaling_factor_3X3);
        h_dummy->GetXaxis()->SetLabelSize(Label_size_3X3*scaling_factor_3X3);
        if(i_pad > 6) h_dummy->GetXaxis()->SetLabelSize(0.075*scaling_factor_3X3);
        h_dummy->GetYaxis()->SetLabelSize(Label_size_3X3*scaling_factor_3X3);
        h_dummy->GetXaxis()->SetTitleSize(Label_size_3X3*scaling_factor_3X3);
        h_dummy->GetYaxis()->SetTitleSize(Label_size_3X3*scaling_factor_3X3);
        h_dummy->GetXaxis()->SetNdivisions(505,'N');
        h_dummy->GetYaxis()->SetNdivisions(505,'N');
        h_dummy->SetLineColor(10);
        h_dummy ->DrawCopy();

        PlotLine(0.0,5,0,0,1,1,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle
        if(i_pad == 7)
        {
            plotTopLegend((char*)"p_{T} (GeV/c)",0.3,0.1,0.1,1,0.0,42,1);
        }
        if(i_pad == 3)
        {
            plotTopLegend((char*)"#it{v}_{2}",0.12,0.5,0.17,1,90,42,1);
        }
    }

    printf("Do_Plot \n");
    DoFillGraph();  // Fill vector vec_tgae with graphs before plotting

    for(Int_t i_entry = 0; i_entry <  N_used_entries; i_entry++)
    {
        Int_t i_select_table = vec_table[i_entry]->GetNumber();
        if(i_select_table == 0) return 0;                           // Select table before plotting
        c_3X3->cd(i_entry+1);
        vec_tgae[i_entry]->SetMarkerSize(0.75);
        vec_tgae[i_entry]->SetMarkerStyle(24);
        vec_tgae[i_entry]->SetMarkerColor(kBlack);
        vec_tgae[i_entry]->Draw("same PZ");
        vec_tgae[i_entry]->Clear();
    }

    return 1;
}
//-------------------------------------------------------------------------------------------------------------------------------------------



//-------------------------------------------------------------------------------------------------------------------------------------------
TestFileList::TestFileList(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h)
{
   //Constructor
   //Create transient frame containing a filelist widget.
   TGLayoutHints *lo;

   fMain = new TGTransientFrame(p, main, w, h);
   fMain->Connect("CloseWindow()", "TestFileList", this, "CloseWindow()");
   fMain->DontCallClose(); // to avoid double deletions.

   //Use hierarchical cleaning
   fMain->SetCleanup(kDeepCleanup);

   TGMenuBar* mb = new TGMenuBar(fMain);
   lo = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1);
   fMain->AddFrame(mb, lo);

   foptions = new TGCompositeFrame(fMain, 60,20,kHorizontalFrame);
   fMain->AddFrame(foptions, lo);

   ffile = new TGCompositeFrame(foptions, 60,20,kVerticalFrame);
   foptions->AddFrame(ffile, lo);

   hep = new TGRadioButton(ffile, "hep file");
   ffile->AddFrame(hep);
   hep->Connect("Clicked()", "TestFileList", this, "hep_clicked()");
   text = new TGRadioButton(ffile, "text file");
   ffile->AddFrame(text);
   text->Connect("Clicked()", "TestFileList", this, "text_clicked()");

   fpt = new TGCompositeFrame(foptions, 60,20,kVerticalFrame);
   foptions->AddFrame(fpt,lo);

   pt  = new TGRadioButton(fpt, "pT");
   fpt->AddFrame(pt);
   pt->Connect("Clicked()", "TestFileList", this, "pt_clicked()");

   pt_high_low  = new TGRadioButton(fpt, "pT_low_high");
   fpt->AddFrame(pt_high_low);
   pt_high_low->Connect("Clicked()", "TestFileList", this, "pt_low_high_clicked()");

   mt  = new TGRadioButton(fpt, "mT");
   fpt->AddFrame(mt);
   mt->Connect("Clicked()", "TestFileList", this, "mt_clicked()");

   mt_low_high  = new TGRadioButton(fpt, "mT_low_high");
   fpt->AddFrame(mt_low_high);
   mt_low_high->Connect("Clicked()", "TestFileList", this, "mt_low_high_clicked()");

   pt_error = new TGCheckButton(fpt, "pT_error");
   fpt->AddFrame(pt_error);

   ferrors = new TGCompositeFrame(foptions, 60,20,kVerticalFrame);
   foptions->AddFrame(ferrors, lo);

   stat_error = new TGCheckButton(ferrors, "stat_error");
   syst_error = new TGCheckButton(ferrors, "syst_error");
   ferrors->AddFrame(stat_error);
   ferrors->AddFrame(syst_error);
   //stat_error->Connect("Clicked()", "TestFileList", this, "stat_error_clicked()");
   syst_error->Connect("Clicked()", "TestFileList", this, "syst_error_clicked()");
   //syst_percent = new TGCheckButton(ferrors,"syst_in_%");
   //ferrors->AddFrame(syst_percent);

   //Text entry for input file name, filled by double click on root files
   Add_file_name = new TGTextEntry(fMain, new TGTextBuffer(15));
   fMain->AddFrame(Add_file_name, lo);

   //When using Open the displayed root file is loaded
   Open = new TGTextButton(fMain, "Open");
   Open->Connect("Clicked()", "TestFileList", this, "DoOpen()");
   fMain->AddFrame(Open, lo);

   //Popup menu to switch view mode
   fMenu = mb->AddPopup("&View");
   fMenu->AddEntry("Lar&ge Icons",kLVLargeIcons);
   fMenu->AddEntry("S&mall Icons",kLVSmallIcons);
   fMenu->AddEntry("&List",       kLVList);
   fMenu->AddEntry("&Details",    kLVDetails);
   fMenu->AddSeparator();
   fMenu->AddEntry("&Close",      10);
   fMenu->Connect("Activated(Int_t)","TestFileList",this,"DoMenu(Int_t)");

   TGListView* lv = new TGListView(fMain, w, h);
   lo = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
   fMain->AddFrame(lv,lo);

   Pixel_t white;
   gClient->GetColorByName("white", white);
   fContents = new TGFileContainer(lv, kSunkenFrame,white);
   fContents->Connect("DoubleClicked(TGFrame*,Int_t)", "TestFileList", this,
                      "OnDoubleClick(TGLVEntry*,Int_t)");

   fMain->CenterOnParent();         // position relative to the parent's window
   fMain->SetWindowName("Data");
   fMain->MapSubwindows();
   fMain->MapWindow();
   fContents->SetDefaultHeaders();
   fContents->DisplayDirectory();
   fContents->AddFile("..");        // up level directory
   fContents->Resize();
   fContents->StopRefreshTimer();   // stop refreshing
   fMain->Resize();
}
//-------------------------------------------------------------------------------------------------------------------------------------------
TestFileList::~TestFileList()
{
   // Cleanup.

   delete fContents;
   fMain->DeleteWindow();  // deletes fMain
}
//-------------------------------------------------------------------------------------------------------------------------------------------



//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::DoOpen()
{
    //Handle Open button
    //Loading input file displayed in text entry

    if (hep->IsDown()) {

        TString file_name = Add_file_name->GetDisplayText();
        newfile = new TFile(file_name.Data());
        printf("open hep \n");
        //delete this;
    }

    if (text->IsDown()) {

        TString file_name = Add_file_name->GetDisplayText();
        string full_name = file_name.Data();
        char* cstr = new char[full_name.length()+1];
        std::strcpy(cstr,full_name.c_str());
        myfile.open((char*)cstr);
        char NoP[50];
        printf("open text \n");
        //delete this;
    }

    flag_file_opened = 1;
}

//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::DoMenu(Int_t mode)
{
   // Switch view mode.

   if (mode<10) {
      fContents->SetViewMode((EListViewMode)mode);
   } else {
      delete this;
   }
}



//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::DisplayFile(const TString &fname)
{
   // Display root files

   Add_file_name->Clear();
   Add_file_name->AppendText(fname); //append file name to text entry
   fContents->SetPagePosition(0,0);
   fContents->SetColHeaders("Name","Title");
   fMain->Resize();
}



//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::DisplayDirectory(const TString &fname)
{
   // Display content of directory.

   fContents->SetDefaultHeaders();
   gSystem->ChangeDirectory(fname);
   fContents->ChangeDirectory(fname);
   fContents->DisplayDirectory();
   fContents->AddFile("..");  // up level directory
   fMain->Resize();
}



//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::DisplayObject(const TString& fname,const TString& name)
{
   // Browse object located in file.


   TDirectory *sav = gDirectory;

   static TFile *file = 0;
   if (file) delete file;     // close
   file = new TFile(fname);   // reopen

   TObject* obj = file->Get(name);
   if (obj) {
       if (!obj->IsFolder())
       {
         obj->Browse(0);
       }
       else obj->Print();
   }
   gDirectory = sav;
}



//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::OnDoubleClick(TGLVEntry *f, Int_t btn)
{
   // Handle double click.

   if (btn != kButton1) return;

   // set kWatch cursor
   ULong_t cur = gVirtualX->CreateCursor(kWatch);
   gVirtualX->SetCursor(fContents->GetId(), cur);

   TString name(f->GetTitle());
   fname = (const char*)f->GetUserData();
   if (fname)
   {
      DisplayObject(fname, name);
   }
   else if (name.EndsWith(".root"))
   {
       DisplayFile(name);
   }
   else if (name.EndsWith(".txt"))
   {
       DisplayFile(name);
   }
   else
   {
      DisplayDirectory(name);
   }
   // set kPointer cursor
   cur = gVirtualX->CreateCursor(kPointer);
   gVirtualX->SetCursor(fContents->GetId(), cur);
}



//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::CloseWindow()
{
   delete this;
}
//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::hep_clicked()
{
    // Handle hep file button
    // See v2_clicked()

    printf("hep_clicked \n");
    text ->SetState(kButtonUp);
}
//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::text_clicked()
{
    // Handle text file button
    // See v2_clicked()

    printf("text_clicked \n");
    hep->SetState(kButtonUp);
}
//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::pt_clicked()
{
    // Handle hep file button
    // See v2_clicked()

    printf("pt_clicked \n");
    pt_high_low ->SetState(kButtonUp);
    mt_low_high ->SetState(kButtonUp);
    mt->SetState(kButtonUp);
}
//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::pt_low_high_clicked()
{
    // Handle text file button
    // See v2_clicked()

    printf("pt_low_high_clicked \n");
    pt->SetState(kButtonUp);
    mt_low_high ->SetState(kButtonUp);
    mt->SetState(kButtonUp);
}
//-------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::mt_clicked()
{
    // Handle hep file button
    // See v2_clicked()

    printf("mt_clicked \n");
    mt_low_high ->SetState(kButtonUp);
    pt->SetState(kButtonUp);
    pt_high_low ->SetState(kButtonUp);
}
//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::mt_low_high_clicked()
{
    // Handle text file button
    // See v2_clicked()

    printf("mt_low_high_clicked \n");
    mt->SetState(kButtonUp);
    pt->SetState(kButtonUp);
    pt_high_low ->SetState(kButtonUp);
}
//-------------------------------------------------------------------------------------------------------------------------------------------
void TestFileList::syst_error_clicked()
{

    printf("syst_error_clicked \n");
    syst_percent = new TGCheckButton(ferrors,"syst_in_%");
    syst_low_high = new TGCheckButton(ferrors,"syst_low_high");
    ferrors->AddFrame(syst_percent);
    ferrors->AddFrame(syst_low_high);
    ferrors->MapSubwindows();
    //ferrors->MapWindow();
    //ferrors->Resize();
    foptions->Resize(400,100);
    fMain->Resize();

}
//-------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------
void SelectData()
{
    new MyMainFrame(gClient->GetRoot(), 20, 40);
}
//--------------------------------------------------------------------------------------------------------------------------------------------


