
#include "./functions_BW.h"


Int_t v2_boost(Double_t Temp_set = 0.1, Double_t rho_0_set = 0.95, Double_t rho_a_set = 0.08,
               Double_t R_x_set = 0.8, Double_t fboost_set = 0.1,
               Int_t flag_loop = 0, Int_t i_rho_a_not_used = 0, Int_t i_R_x = 0, Int_t i_fboost = 0)
{

    // .x v2_boost.cc++(0.15,0.95,0.08,0.8,0.1,0,0,0,0)
    // .x v2_boost.cc++(0.14,0.925,0.15,0.9,1.0,0,0,0,0)

    //--------------------------------------------------------------------
    cout << "v2_boost started" << endl;

    TGaxis::SetMaxDigits(3);
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(0);
    //--------------------------------------------------------------------

    init_data();
    Init_density();
    Init_v2_Mathematica();
    init_pT_spectra_data();
    make_5_60_spectra();

    //--------------------------------------------------------------------
    TFile* outputfile = new TFile(Form("out_v2_boost_R_x%d_fb%d.root",i_R_x,i_fboost),"RECREATE");
    //--------------------------------------------------------------------


    //--------------------------------------------------------------------
    TF1* f_FlowFitFunc = new TF1("f_FlowFitFunc",FlowFitFunc,0,10,5);
    //--------------------------------------------------------------------



    //--------------------------------------------------------------------
    const Int_t    N_quarks         = 100; // 2000
    const Double_t x_fac            = 0.7;  // 1.0
    const Double_t y_fac            = 1.0;  // 2.5
    const Double_t R_scale          = R_Pb*0.75;  // 1.0
    const Double_t beta_surface     = 0.75;  // 0.7
    const Double_t temperature      = 0.075; // 0.02 0.25, 0.07
    const Int_t    N_events         = 50000; // 50000
    const Double_t rel_distance_max = 0.2; // 0.5
    const Double_t rel_momentum_max = 0.5; // 0.3

    Double_t rho_0 = 0.89;
    Double_t rho_a = 0.05; // 3.02

    //--------------------------------------------------------------------



    //--------------------------------------------------------------------
    TPolyLine* pl_geometric_shape = new TPolyLine();
    Double_t x_pos       = 0.0;
    Double_t y_pos       = 0.0;
    Double_t x_pos_point = 0.0;
    Double_t y_pos_point = 0.0;
    Double_t z_pos_point = 0.0;
    for(Double_t phi = 0.0; phi <= 360.0; phi+= 1.0)
    {
        Double_t radius = geometric_shape(R_scale,phi,x_fac,y_fac,x_pos,y_pos);
        pl_geometric_shape ->SetNextPoint(x_pos,y_pos);
    }

    Double_t size_hist = 1.5*R_scale*TMath::Max(x_fac,y_fac);
    h2D_geometric_shape = new TH2D("h2D_geometric_shape","h2D_geometric_shape",100,-size_hist,size_hist,100,-size_hist,size_hist);

    get_geometric_shape(h2D_geometric_shape,x_fac,y_fac,R_scale);
    //--------------------------------------------------------------------



    //--------------------------------------------------------------------
    //TF1* f_FlowFitFunc = new TF1("f_FlowFitFunc",PtFitFunc,0,10,1);
    f_LevyFitFunc         = new TF1("f_LevyFitFunc",PtFitFunc2_mod,0.0,10.0,4);
    f_LevyFitFunc ->SetParameter(1,temperature); // temperature
    f_LevyFitFunc ->SetParameter(2,1.0);
    f_LevyFitFunc ->SetParameter(3,0.0);
    //f_LevyFitFunc ->Draw();

    f_JetPtFunc         = new TF1("f_JetPtFunc",PtFitFunc2_mod,0.0,10.0,4);
    f_JetPtFunc ->SetParameter(0,0.2); // mass
    f_JetPtFunc ->SetParameter(1,0.9); // temperature 2.9
    f_JetPtFunc ->SetParameter(2,1.0);
    f_JetPtFunc ->SetParameter(3,0.0);

    TF1* f_PtFitFunc         = new TF1("f_PtFitFunc",PtFitFunc,0.0,10.0,1);
    f_PtFitFunc ->SetParameter(0,temperature);


    TCanvas* c_pT_functions = new TCanvas("c_pT_functions","c_pT_functions",100,10,500,500);
    c_pT_functions->cd()->SetRightMargin(0.20);
    c_pT_functions->cd()->SetTopMargin(0.08);
    c_pT_functions->cd();
    for(Int_t i_quark_mass = 0; i_quark_mass < N_masses; i_quark_mass++)
    {
        Double_t quark_mass = arr_quark_mass_meson[i_quark_mass];
        f_LevyFitFunc ->SetParameter(0,quark_mass);
        f_LevyFitFunc ->SetLineColor(arr_color_mass[i_quark_mass]);
        if(i_quark_mass == 0) f_LevyFitFunc ->DrawCopy();
        else f_LevyFitFunc ->DrawCopy("same");
    }
    //--------------------------------------------------------------------



    //--------------------------------------------------------------------
    TCanvas* c_geometric_shape = Draw_2D_histo_and_canvas(h2D_geometric_shape,"c_geometric_shape",1010,820,0.0,-1.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    c_geometric_shape->cd()->SetRightMargin(0.20);
    c_geometric_shape->cd()->SetTopMargin(0.08);
    c_geometric_shape->cd();
    pl_geometric_shape ->SetLineColor(kRed);
    pl_geometric_shape ->SetLineWidth(2);
    pl_geometric_shape ->SetLineStyle(1);
    pl_geometric_shape ->Draw();

    Double_t x_pos_nucl_A = R_Pb*0.3;
    Draw_Circle_2D_new(R_Pb,R_Pb+0.1,1,20,kRed,1,1,x_pos_nucl_A,0.0,kRed,0.5,0.5);
    Draw_Circle_2D_new(R_Pb,R_Pb+0.1,1,20,kBlue,1,1,-x_pos_nucl_A,0.0,kBlue,0.5,0.5);

    PlotLine(-2.5,2.5,0.0,0.0,kBlack,1,1); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
    PlotLine(0.0,0.0,-2.8,2.8,kBlack,1,1); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)



    TPolyMarker* quark_pm = new TPolyMarker();
    vector<TArrow*> vec_beta;
    vector<TArrow*> vec_beta_perp;
    Int_t size = (Int_t)vec_beta.size();

    vector< vector<TVector3> > vec_vec_surf_perp = geometric_shape_perp(R_scale*1.0,x_fac,y_fac);
    vector<TArrow*> vec_arrow_surf_perp;
    for(Int_t i_perp = 0; i_perp < (Int_t)vec_vec_surf_perp[0].size(); i_perp++)
    {
        vec_arrow_surf_perp.push_back(new TArrow(vec_vec_surf_perp[0][i_perp].X(),vec_vec_surf_perp[0][i_perp].Y(),
                                                 vec_vec_surf_perp[0][i_perp].X()+vec_vec_surf_perp[1][i_perp].X(),
                                                 vec_vec_surf_perp[0][i_perp].Y()+vec_vec_surf_perp[1][i_perp].Y(),0.01,"|>"));
    }

    for(Int_t i_perp = 0; i_perp < (Int_t)vec_arrow_surf_perp.size(); i_perp++)
    {
        vec_arrow_surf_perp[i_perp]->SetAngle(30.0);
        vec_arrow_surf_perp[i_perp]->SetLineWidth(2);
        vec_arrow_surf_perp[i_perp]->SetLineColor(kRed);
        vec_arrow_surf_perp[i_perp]->SetFillColor(kRed);
        //vec_arrow_surf_perp[i_perp]->Draw();
    }

    for(Int_t i_quark = 0; i_quark < 50; i_quark++)
    {
        h2D_geometric_shape ->GetRandom2(x_pos_point,y_pos_point);
        quark_pm ->SetNextPoint(x_pos_point,y_pos_point);

        Double_t radius_point = TMath::Sqrt(x_pos_point*x_pos_point + y_pos_point*y_pos_point);
        Double_t phi          = TMath::RadToDeg()*TMath::ATan2(y_pos_point/(y_fac*R_scale),x_pos_point/(x_fac*R_scale));
        TVector3 tv3_boost      = get_boost_vector(R_scale,x_fac,y_fac,x_pos_point,y_pos_point,0.0,rho_0,rho_a);
        Double_t quark_pT          = f_LevyFitFunc->GetRandom();
        Double_t quark_thermal_phi = ran.Rndm()*360.0;
        Double_t quark_thermal_eta = (ran.Rndm()-0.5)*2.0;
        quark_thermal_eta = 0.0;
        TLorentzVector tlv_quark;
        tlv_quark.SetPtEtaPhiM(quark_pT,quark_thermal_eta,quark_thermal_phi*TMath::DegToRad(),0.14); // thermal

        tlv_quark.Boost(tv3_boost);

        vec_beta.push_back(new TArrow(x_pos_point,y_pos_point,x_pos_point+tv3_boost.X(),y_pos_point+tv3_boost.Y(),0.01,"|>"));
    }


    for(Int_t i_quark = 0; i_quark < (Int_t)vec_beta.size(); i_quark++)
    {
        vec_beta[i_quark]->SetAngle(30.0);
        vec_beta[i_quark]->SetLineWidth(2);
        vec_beta[i_quark]->SetLineColor(kRed);
        vec_beta[i_quark]->SetFillColor(kRed);
        vec_beta[i_quark]->Draw();
    }

    /*
    for(Int_t i_quark = 0; i_quark < (Int_t)vec_beta_perp.size(); i_quark++)
    {
        vec_beta_perp[i_quark]->SetAngle(30.0);
        vec_beta_perp[i_quark]->SetLineWidth(2);
        vec_beta_perp[i_quark]->SetLineColor(kRed);
        vec_beta_perp[i_quark]->SetFillColor(kRed);
        vec_beta_perp[i_quark]->Draw();
    }
    */

    quark_pm ->SetMarkerColor(kBlue);
    quark_pm ->SetMarkerSize(0.2);
    quark_pm ->SetMarkerStyle(20);
    quark_pm ->Draw();
    //--------------------------------------------------------------------



    //--------------------------------------------------------------------
    h2D_density_Glauber ->GetXaxis()->SetTitle("x (fm)");
    h2D_density_Glauber ->GetYaxis()->SetTitle("y (fm)");
    TCanvas* c_density_Glauber = Draw_2D_histo_and_canvas(h2D_density_Glauber,"c_density_Glauber",1010,820,0.0,-1.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    c_density_Glauber->cd()->SetRightMargin(0.20);
    c_density_Glauber->cd()->SetTopMargin(0.08);
    c_density_Glauber->cd();

    pl_geometric_shape ->Draw();
    Draw_Circle_2D_new(R_Pb,R_Pb+0.1,1,20,kRed,1,1,x_pos_nucl_A,0.0,kRed,1.0,0.1);
    Draw_Circle_2D_new(R_Pb,R_Pb+0.1,1,20,kGreen,1,1,-x_pos_nucl_A,0.0,kGreen,1.0,0.1);
    //--------------------------------------------------------------------



    //--------------------------------------------------------------------
    vector<TProfile*> tp_v2_vs_pT_mesons;
    vector<TH1F*> h_dN_dpT_jets;
    tp_v2_vs_pT_mesons.resize(N_masses);
    h_dN_dpT_mesons.resize(N_masses);
    h_dN_dpT_jets.resize(N_masses);

    // pi, K, p, phi, Omega, D0, J/Psi, Upsilon
    Double_t high_pt_range[8]         = {5.0,5.0,6.0,6.0,7.0,10.0,15.0,15.0};
    Int_t    N_bins_v2_pt_range[8]    = {25,25,25,25,25,25,25,25};
    Int_t    N_bins_dNdpT_pt_range[8] = {100,100,100,100,100,100,100,100};

    for(Int_t i_mass = 0; i_mass < N_masses; i_mass++)
    {
        tp_v2_vs_pT_mesons[i_mass]  = new TProfile(Form("tp_v2_vs_pT_mesons_%d",i_mass),Form("tp_v2_vs_pT_mesons_%d",i_mass),N_bins_v2_pt_range[i_mass],0,high_pt_range[i_mass]);
        h_dN_dpT_mesons[i_mass]     = new TH1F(Form("h_dN_dpT_mesons_%d",i_mass),Form("h_dN_dpT_mesons_%d",i_mass),N_bins_dNdpT_pt_range[i_mass],0,high_pt_range[i_mass]);
        h_dN_dpT_jets[i_mass]       = new TH1F(Form("h_dN_dpT_jets_%d",i_mass),Form("h_dN_dpT_jets_%d",i_mass),N_bins_dNdpT_pt_range[i_mass],0,high_pt_range[i_mass]);
    }
    TProfile* tp_v2_vs_pT_thermic = new TProfile("tp_v2_vs_pT_thermic","tp_v2_vs_pT_thermic",150,0,3.0);

    TH1F* tp_dN_dphi_vs_phi   = new TH1F("tp_dN_dphi_vs_phi","tp_dN_dphi_vs_phi",100,-TMath::Pi(),TMath::Pi());
    TH2D* h2D_pos_XY_neg_v2   = new TH2D("h2D_pos_XY_neg_v2","h2D_pos_XY_neg_v2",100,-3,3,100,-3,3);
    TH1F* h_pos_XY_neg_v2_pT  = new TH1F("h_pos_XY_neg_v2_pT","h_pos_XY_neg_v2_pT",100,0,10.0);
    TH2D* h2D_pos_XY_pos_v2   = new TH2D("h2D_pos_XY_pos_v2","h2D_pos_XY_pos_v2",100,-3,3,100,-3,3);
    TH1F* h_pos_XY_pos_v2_pT  = new TH1F("h_pos_XY_pos_v2_pT","h_pos_XY_pos_v2_pT",100,0,10.0);
    TH2D* h2D_rapidity_vs_eta = new TH2D("h2D_rapidity_vs_eta","h2D_rapidity_vs_eta",100,-2.0,2.0,100,-2.0,2.0);
    TH2D* h2D_cos_theta_vs_pT = new TH2D("h2D_cos_theta_vs_pT","h2D_cos_theta_vs_pT",100,0,5,100,-1,1);
    //--------------------------------------------------------------------





    //--------------------------------------------------------------------
#if 0
    Double_t rho_0_start = 0.3;  // 0.7
    Double_t rho_0_stop  = 1.4; // 1.2
    Double_t delta_rho_0 = 0.1;  // 0.1

    Double_t rho_a_start = 0.0; // 0.0
    Double_t rho_a_stop  = 0.3; // 0.1
    Double_t delta_rho_a = 0.02; // 0.01

    Double_t R_x_start   = 0.2; // 0.3
    Double_t R_x_stop    = 0.9; // 0.6
    Double_t delta_R_x   = 0.1; // 0.1

    Double_t Temp_start   = 0.05; // 0.05
    Double_t Temp_stop    = 0.2; // 0.2
    Double_t delta_Temp   = 0.025; // 0.025
#endif

#if 1
    Double_t rho_0_start = 0.5;  // 0.7
    Double_t rho_0_stop  = 0.9; // 1.2
    Double_t delta_rho_0 = 0.1;  // 0.1

    Double_t rho_a_start = 0.0; // 0.0
    Double_t rho_a_stop  = 0.3; // 0.1
    Double_t delta_rho_a = 0.02; // 0.01

    Double_t R_x_start   = 0.2; // 0.3
    Double_t R_x_stop    = 0.9; // 0.6
    Double_t delta_R_x   = 0.1; // 0.1

    Double_t Temp_start   = 0.1; // 0.05
    Double_t Temp_stop    = 0.21; // 0.2
    Double_t delta_Temp   = 0.025; // 0.025
#endif

#if 0
    // Good v2 fits for trans -> long boost order
    Double_t Temp_best  = 0.075; // 0.075
    Double_t rho_0_best = 1.3; // 0.9
    Double_t rho_a_best = 0.08; // 0.05
    Double_t R_x_best   = 0.8; // 0.4
#endif

#if 0
    //
    Double_t Temp_best  = 0.175; // 0.075
    Double_t rho_0_best = 0.95; // 0.9
    Double_t rho_a_best = 0.08; // 0.05
    Double_t R_x_best   = 0.8; // 0.4
#endif

#if 0
    // From Markus? Good J/Psi spectra fit
    Double_t Temp_best  = 0.1565; // 0.075
    Double_t rho_0_best = 0.69314718; // 0.9
    Double_t rho_a_best = 0.18; // 0.05
    Double_t R_x_best   = 0.5; // 0.4
#endif

#if 0
    // Good v2 fits for long -> trans boost order
    Double_t Temp_best  = 0.175; // 0.075
    Double_t rho_0_best = 1.2; // 0.9
    Double_t rho_a_best = 0.18; // 0.05
    Double_t R_x_best   = 0.7; // 0.4
#endif

#if 0
    // Good v2 fits for long -> trans boost order
    Double_t Temp_best  = 0.2; // 0.075
    Double_t rho_0_best = 0.8; // 0.9
    Double_t rho_a_best = 0.04; // 0.05
    Double_t R_x_best   = 0.8; // 0.4
#endif


#if 0
    Double_t Temp_best  = 0.165; // 0.2
    Double_t rho_0_best = 0.6; // 0.8
    Double_t rho_a_best = 0.1; // 0.2
    Double_t R_x_best   = 0.7; // 0.8
#endif

    // External settings
#if 1
    Double_t Temp_best   = Temp_set;
    Double_t rho_0_best  = rho_0_set;
    Double_t rho_a_best  = rho_a_set;
    Double_t R_x_best    = R_x_set;
    Double_t fboost_best = fboost_set;
#endif


    // itt: 795, chi2_best: 0.103, rho_0: 1.200, rho_a: 0.000, R_x: 0.600  only pions
    // itt: 823, chi2_best: 0.396, rho_0: 1.300, rho_a: 0.260, R_x: 0.600  pions and kaons
    // itt: 1147, chi2_best: 0.439, rho_0: 1.300, rho_a: 0.140, R_x: 0.800 all
    // itt: 2464, chi2_best: 0.397, Temp: 0.075, rho_0: 1.300, rho_a: 0.080, R_x: 0.800, s2: 0.11
    // itt: 7569, chi2_best: 2.255, Temp: 0.1750, rho_0: 1.200, rho_a: 0.180, R_x: 0.700   all particle, 3 GeV/c
    // itt: 5040, chi2_best: 204.894, Temp: 0.1250, rho_0: 0.900, rho_a: 0.000, R_x: 0.800 -> good fit to pi,K,p,J/Psi spectra, only J/Psi were used

    // itt: 7906, chi2_best: 418.796, Temp: 0.1750, rho_0: 1.300, rho_a: 0.020, R_x: 0.900 -> long, trans, J/Psi spectra and v2
    // itt: 8987, chi2_best: 51.900, Temp: 0.2000, rho_0: 0.800, rho_a: 0.040, R_x: 0.800  -> trans, long, J/Psi spectra and v2

    // itt: 96, chi2_best: 450.598, Temp: 0.2000, rho_0: 0.800, rho_a: 0.120, R_x: 0.800, 0.2 ratio for long -> trans, quite good spectra, v2 not too far off
    // itt: 395, chi2_best: 280.579, Temp: 0.2000, rho_0: 0.600, rho_a: 0.100, R_x: 0.700 same as above






    TProfile* tp_v2_vs_pT_opt_best = NULL;
    h_dNdpT_best = NULL;

#if 0
    printf("Optimize_v2 started");
    Long64_t N_particles = 1000000; // 10000
    Optimize_v2(rho_0_start, rho_0_stop, delta_rho_0,
                rho_a_start, rho_a_stop, delta_rho_a,
                R_x_start, R_x_stop, delta_R_x,
                Temp_start, Temp_stop, delta_Temp,
                N_particles, R_scale,
                rho_0_best, rho_a_best, R_x_best, Temp_best,
                tp_v2_vs_pT_opt_best, h_dNdpT_best
               );
#endif

    //get_geometric_shape(h2D_geometric_shape,R_x_best*0.88,1.0,R_scale);

    Double_t Temp_loop_start  = Temp_best;
    Double_t rho_0_loop_start = rho_0_best;
    Double_t rho_a_loop_start = rho_a_best;
    Int_t N_Temp_loop  = 1;
    Int_t N_rho_0_loop = 1;
    Int_t N_rho_a_loop = 1;

    if(flag_loop)
    {
        Temp_loop_start  = 0.08;
        rho_0_loop_start = 0.3;
        rho_a_loop_start = 0.0;
        N_Temp_loop  = 9; // 9
        N_rho_0_loop = 9; // 9
        N_rho_a_loop = 9; // 9

        //Double_t arr_rho_a[9]   = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4};
        //Double_t arr_R_x[9]     = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
        //Double_t arr_f_boost[9] = {0.0,0.05,0.1,0.15,0.2,0.4,0.6,0.8,1.0};

        //rho_a_best  = arr_rho_a[i_rho_a];
        R_x_best    = arr_R_x[i_R_x];
        fboost_best = arr_f_boost[i_fboost];
    }

    TList* list_profiles = new TList();
    list_profiles ->SetOwner(kTRUE);
    list_profiles ->SetName(Form("list_BW_Rx%d_fb%d",i_R_x,i_fboost));

    for(Int_t i_Temp = 0; i_Temp < N_Temp_loop; i_Temp++)
    {
        printf("i_Temp: %d \n",i_Temp);
        if(flag_loop) Temp_best = Temp_loop_start + i_Temp*0.02;
        for(Int_t i_rho_0 = 0; i_rho_0 < N_rho_0_loop; i_rho_0++)
        {
            printf("i_rho_0: %d \n",i_rho_0);
            if(flag_loop) rho_0_best = rho_0_loop_start + i_rho_0*0.125;

            for(Int_t i_rho_a = 0; i_rho_a < N_rho_a_loop; i_rho_a++)
            {
                printf("i_rho_a: %d \n",i_rho_a);
                if(flag_loop) rho_a_best = rho_a_loop_start + i_rho_a*0.05;

                get_geometric_shape(h2D_geometric_shape,R_x_best,1.0,R_scale);
                f_LevyFitFunc ->SetParameter(1,Temp_best);

                Double_t s2 = get_s2(R_x_best,1.0);
                printf("R_x_best: %4.3f, s2: %4.3f \n",R_x_best,s2);
                //--------------------------------------------------------------------

                printf("Temp \n");
                //--------------------------------------------------------------------
                //for(Int_t i_quark_mass = 0; i_quark_mass < N_masses; i_quark_mass++)
                for(Int_t i_quark_mass = 0; i_quark_mass < 5; i_quark_mass++)
                {
                    Double_t quark_mass = arr_quark_mass_meson[i_quark_mass];
                    printf("--------------- i_quark_mass: %d, mass: %4.3f GeV --------------- \n",i_quark_mass,quark_mass);
                    f_LevyFitFunc ->SetParameter(0,quark_mass);
                    f_JetPtFunc   ->SetParameter(0,quark_mass); // mass
                    Long64_t N_events_use = N_events;
                    if(i_quark_mass > 2) N_events_use = 1*N_events;
                    for(Long64_t i_event = 0; i_event < N_events_use; i_event++)
                    {
                        //--------------------------------------------------------------------
                        // Create single event quark distributions
                        //printf("Create single event quark distributions \n");

                        if (i_event != 0  &&  i_event % 1000 == 0)
                            cout << "." << flush;
                        if (i_event != 0  &&  i_event % 10000 == 0)
                        {
                            cout << 100.0*((Double_t)i_event)/N_events_use << "\% processed" << endl;
                        }

                        for(Int_t i_quark = 0; i_quark < N_quarks; i_quark++)
                        {
                            //h2D_density_Glauber ->GetRandom2(x_pos_point,y_pos_point);
                            h2D_geometric_shape ->GetRandom2(x_pos_point,y_pos_point);  // sample postion of particle
                            Double_t weight = 1.0;
                            //if(fabs(x_pos_point) > 2.4) weight = 0.0;
                            //weight = TMath::Power(TMath::Gaus(x_pos_point,0.0,2.5),4.0);

                            // Sample z-rapidity
                            Double_t z_rapidity = (ran.Rndm()-0.5)*2.0; // get bjorken z-rapidity
                            //Double_t z_rapidity = ran.Gaus(0.0,1.2);

                            Double_t beta_z = (TMath::Exp(2.0*z_rapidity) - 1.0)/(TMath::Exp(2.0*z_rapidity) + 1.0); // = TMath::TanH(z_rapidity), z-beta
                            //printf("beta_z: %4.3f, TanH(y): %4.3f \n,",beta_z,TMath::TanH(z_rapidity));

                            Double_t quark_thermal_phi       = ran.Rndm()*360.0; // sample phi
                            Double_t quark_thermal_cos_theta = (ran.Rndm()-0.5)*2.0; // sample cos theta -> polar angle
                            //Double_t quark_thermal_theta     = TMath::Sign(1.0,ran.Rndm()-0.5)*TMath::ACos(quark_thermal_cos_theta);
                            Double_t quark_thermal_theta     = TMath::ACos(quark_thermal_cos_theta); // get theta
                            Double_t quark_thermal_eta       = -TMath::Log(TMath::Tan(quark_thermal_theta/2.0)); // get eta
                            //Double_t quark_pT  = ran.Rndm()*4.0;
                            //Double_t quark_pT_weight = f_LevyFitFunc->Eval(quark_pT);
                            Double_t quark_pT                = f_LevyFitFunc->GetRandom()*TMath::Sin(quark_thermal_theta); // sample p, transform to pT
                            //Double_t quark_pT                = f_LevyFitFunc->GetRandom();
                            //Double_t quark_pT                = f_PtFitFunc->GetRandom();


                            TLorentzVector tlv_quark;
                            tlv_quark.SetPtEtaPhiM(quark_pT,quark_thermal_eta,quark_thermal_phi*TMath::DegToRad(),quark_mass); // thermal

                            //if(quark_thermal_theta < 0.1)
                            //{
                            //    printf("phi: %4.3f, theta: %4.3f, eta: %4.3f, pT: %4.3f \n",quark_thermal_phi,quark_thermal_theta,quark_thermal_eta,quark_pT);
                            //    tlv_quark.Print();
                            //}


                            Double_t pT_thermic = tlv_quark.Pt();

                            TVector3  tv3_boost  = get_boost_vector(R_scale,R_x_best,1.0,x_pos_point,y_pos_point,z_rapidity,rho_0_best,rho_a_best); // get transverse boost vector
                            //tv3_boost.SetZ(beta_z);
                            //tv3_boost.SetZ(0.0);

#if 0
                            tv3_boost.SetX(-tv3_boost.X());
                            tv3_boost.SetY(-tv3_boost.Y());
                            tv3_boost.SetZ(-tv3_boost.Z());
                            TVector3  tv3_boost_long(0.0,0.0,-beta_z); // longitudinal boost vector
#endif

                            TVector3  tv3_boost_long(0.0,0.0,beta_z); // longitudinal boost vector

                            // boost order is important! -> work in progress


                             TMatrixD TL_Matrix(4,4);
                             TMatrixD TL_MatrixInv(4,4);
                             TMatrixD TL_MatrixBInv(4,4);
                             TMatrixD TL_MatrixMult(4,4);

                             Double_t Yb    = TMath::ATanH(beta_z);
                             Double_t rhob  = TMath::ATanH(tv3_boost.Mag());
                             Double_t phi_s = TMath::ATan2(y_pos_point,x_pos_point);;
                             Double_t phi_b = TMath::ATan(TMath::Tan(phi_s)/TMath::Power(1.0/R_x_best,2.0));

                             //----------------------------
                             // transformation from lab to cell as used by Klaus
                             TL_Matrix[0][0] = TMath::CosH(Yb)*TMath::CosH(rhob);
                             TL_Matrix[0][1] = -TMath::Cos(phi_b)*TMath::SinH(rhob);
                             TL_Matrix[0][2] = -TMath::Sin(phi_b)*TMath::SinH(rhob);
                             TL_Matrix[0][3] = -TMath::CosH(rhob)*TMath::SinH(Yb);

                             TL_Matrix[1][0] = -TMath::CosH(Yb)*TMath::SinH(rhob);
                             TL_Matrix[1][1] = TMath::Cos(phi_b)*TMath::CosH(rhob);
                             TL_Matrix[1][2] = TMath::CosH(rhob)*TMath::Sin(phi_b);
                             TL_Matrix[1][3] = TMath::SinH(Yb)*TMath::SinH(rhob);

                             TL_Matrix[2][0] = 0.0;
                             TL_Matrix[2][1] = -TMath::Sin(phi_b);
                             TL_Matrix[2][2] = TMath::Cos(phi_b);
                             TL_Matrix[2][3] = 0.0;

                             TL_Matrix[3][0] = -TMath::SinH(Yb);
                             TL_Matrix[3][1] = 0.0;
                             TL_Matrix[3][2] = 0.0;
                             TL_Matrix[3][3] = TMath::CosH(Yb);
                             //----------------------------


                             //----------------------------
                             // Inverse transformation is identical to my trans. first then long. transformation
                             TL_MatrixInv[0][0] = TMath::CosH(Yb)*TMath::CosH(rhob);
                             TL_MatrixInv[0][1] = TMath::CosH(Yb)*TMath::SinH(rhob);
                             TL_MatrixInv[0][2] = 0.0;
                             TL_MatrixInv[0][3] = TMath::SinH(Yb);

                             TL_MatrixInv[1][0] = TMath::Cos(phi_b)*TMath::SinH(rhob);
                             TL_MatrixInv[1][1] = TMath::Cos(phi_b)*TMath::CosH(rhob);
                             TL_MatrixInv[1][2] = -TMath::Sin(phi_b);
                             TL_MatrixInv[1][3] = 0.0;

                             TL_MatrixInv[2][0] = TMath::Sin(phi_b)*TMath::SinH(rhob);
                             TL_MatrixInv[2][1] = TMath::CosH(rhob)*TMath::Sin(phi_b);
                             TL_MatrixInv[2][2] = TMath::Cos(phi_b);
                             TL_MatrixInv[2][3] = 0.0;

                             TL_MatrixInv[3][0] = TMath::CosH(rhob)*TMath::SinH(Yb);
                             TL_MatrixInv[3][1] = TMath::SinH(Yb)*TMath::SinH(rhob);
                             TL_MatrixInv[3][2] = 0.0;
                             TL_MatrixInv[3][3] = TMath::CosH(Yb);
                             //----------------------------


                            //----------------------------
                            /*
                             // Inverse transformation for radial and then long in Klaus's case,
                             TL_MatrixBInv[0][0] = TMath::CosH(Yb)*TMath::CosH(rhob);
                             TL_MatrixBInv[0][1] = -TMath::Cos(phi_b)*TMath::CosH(Yb)*TMath::SinH(rhob);
                             TL_MatrixBInv[0][2] = -TMath::CosH(Yb)*TMath::Sin(phi_b)*TMath::SinH(rhob);
                             TL_MatrixBInv[0][3] = -TMath::SinH(Yb);

                             TL_MatrixBInv[1][0] = -TMath::SinH(rhob);
                             TL_MatrixBInv[1][1] = TMath::Cos(phi_b)*TMath::CosH(rhob);
                             TL_MatrixBInv[1][2] = TMath::CosH(rhob)*TMath::Sin(phi_b);
                             TL_MatrixBInv[1][3] = 0.0;

                             TL_MatrixBInv[2][0] = 0.0;
                             TL_MatrixBInv[2][1] = -TMath::Sin(phi_b);
                             TL_MatrixBInv[2][2] = TMath::Cos(phi_b);
                             TL_MatrixBInv[2][3] = 0.0;

                             TL_MatrixBInv[3][0] = -TMath::CosH(rhob)*TMath::SinH(Yb);
                             TL_MatrixBInv[3][1] = TMath::Cos(phi_b)*TMath::SinH(Yb)*TMath::SinH(rhob);
                             TL_MatrixBInv[3][2] = TMath::Sin(phi_b)*TMath::SinH(Yb)*TMath::SinH(rhob);
                             TL_MatrixBInv[3][3] = TMath::CosH(Yb);
                             */

                            /*
                             // Inverse transformation for radial and then long in Klaus's case,
                             TL_MatrixBInv[0][0] = TMath::CosH(Yb)*TMath::CosH(rhob);
                             TL_MatrixBInv[0][1] = TMath::SinH(rhob);
                             TL_MatrixBInv[0][2] = 0.0;
                             TL_MatrixBInv[0][3] = TMath::CosH(rhob)*TMath::SinH(Yb);

                             TL_MatrixBInv[1][0] = TMath::Cos(phi_b)*TMath::CosH(Yb)*TMath::SinH(rhob);
                             TL_MatrixBInv[1][1] = TMath::Cos(phi_b)*TMath::CosH(rhob);
                             TL_MatrixBInv[1][2] = -TMath::Sin(phi_b);
                             TL_MatrixBInv[1][3] = TMath::Cos(phi_b)*TMath::SinH(Yb)*TMath::SinH(rhob);

                             TL_MatrixBInv[2][0] = TMath::CosH(Yb)*TMath::Sin(phi_b)*TMath::SinH(rhob);
                             TL_MatrixBInv[2][1] = TMath::CosH(rhob)*TMath::Sin(phi_b);
                             TL_MatrixBInv[2][2] = TMath::Cos(phi_b);
                             TL_MatrixBInv[2][3] = TMath::Sin(phi_b)*TMath::SinH(Yb)*TMath::SinH(rhob);

                             TL_MatrixBInv[3][0] = TMath::SinH(Yb);
                             TL_MatrixBInv[3][1] = 0.0;
                             TL_MatrixBInv[3][2] = 0.0;
                             TL_MatrixBInv[3][3] = TMath::CosH(Yb);
                             //----------------------------

                             TL_MatrixMult = TL_MatrixInv*TL_Matrix;
                             //TL_MatrixMult.Print();


                             //TVectorD vec_Lprime = TL_Matrix*vec_L;
                             //TVectorD vec_Lprime = TL_MatrixInv*vec_L;
                             */

                             TVectorD vec_L(4);
                             vec_L[0] = tlv_quark.E();
                             vec_L[1] = tlv_quark.Px();
                             vec_L[2] = tlv_quark.Py();
                             vec_L[3] = tlv_quark.Pz();
                             //TVectorD vec_Lprime = TL_MatrixInv*vec_L;


                             Init_Lorentz_Matrices(Yb,phi_b,rhob);
#if 0
                             //-------------------------------------------------------
                             // LR^-1TR * vec
                             TVectorD vec_Lprime = TL_Matrix_Lambda_R*vec_L;
                             vec_L = vec_Lprime;
                             vec_Lprime = TL_Matrix_Lambda_T*vec_L;
                             vec_L = vec_Lprime;
                             vec_Lprime = TL_Matrix_Lambda_RInv*vec_L;
                             vec_L = vec_Lprime;
                             vec_Lprime = TL_Matrix_Lambda_L*vec_L;
                             //-------------------------------------------------------
#endif


#if 1
                             //-------------------------------------------------------
                             // R^-1TLR * vec
                             TVectorD vec_Lprime = TL_Matrix_Lambda_R*vec_L;
                             vec_L = vec_Lprime;
                             vec_Lprime = TL_Matrix_Lambda_L*vec_L;
                             vec_L = vec_Lprime;
                             vec_Lprime = TL_Matrix_Lambda_T*vec_L;
                             vec_L = vec_Lprime;
                             vec_Lprime = TL_Matrix_Lambda_RInv*vec_L;
                             //-------------------------------------------------------
#endif

                             //tlv_quark.SetPxPyPzE(vec_Lprime[1],vec_Lprime[2],vec_Lprime[3],vec_Lprime[0]);



#if 0
                            tlv_quark.Boost(tv3_boost);
                            tlv_quark.Boost(tv3_boost_long);
                            //tlv_quark.Boost(tv3_boost);
#endif

#if 1
                            if(ran.Rndm() > fboost_best)
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

                            if(i_quark_mass == 0) h2D_rapidity_vs_eta ->Fill(eta,rapidity);

                            /*
                             // Standard acceptance cuts
                             if(i_quark_mass < 3  && fabs(eta) > 1.0) continue;
                             if(i_quark_mass < 3  && fabs(rapidity) > 0.5) continue;
                             if(i_quark_mass >= 3 && (rapidity < 2.5 || rapidity > 4.0)) continue;
                             */

                            //if(i_quark_mass < 5  && fabs(eta) > 0.1) continue;
                            //if(i_quark_mass < 5  && fabs(rapidity) > 0.5) continue;
                            if(i_quark_mass < 5  && fabs(rapidity) > 0.1) continue;
                            //if(i_quark_mass >= 3 && (rapidity < 2.5 || rapidity > 4.0)) continue;

                            h2D_cos_theta_vs_pT ->Fill(quark_pT,quark_thermal_cos_theta);

                            // Sample jet particles
                            Double_t jet_pT                = f_JetPtFunc->GetRandom()*TMath::Sin(quark_thermal_theta); // sample p, transform to pT
                            TLorentzVector tlv_jet;
                            tlv_jet.SetPtEtaPhiM(jet_pT,quark_thermal_eta,quark_thermal_phi*TMath::DegToRad(),quark_mass); // jet particle
                            Double_t cos_phin_jet = TMath::Cos(2.0*tlv_jet.Phi());
                            Double_t weight_jet = 0.0;

                            if(i_quark_mass == 2 && pT_lab > 3.5 && pT_lab < 4.5) tp_dN_dphi_vs_phi ->Fill(tlv_quark.Phi());
                            //printf("phi: %4.3f \n",tlv_quark.Phi());

                            tp_v2_vs_pT_mesons[i_quark_mass] ->Fill(pT_lab,cos_phin,weight);
                            tp_v2_vs_pT_mesons[i_quark_mass] ->Fill(jet_pT,cos_phin_jet,weight_jet);

                            h_dN_dpT_mesons[i_quark_mass]    ->Fill(pT_lab,weight);
                            h_dN_dpT_jets[i_quark_mass]      ->Fill(jet_pT,weight_jet);
                            if(i_quark_mass == 2 && cos_phin < -0.4)
                            {
                                h2D_pos_XY_neg_v2  ->Fill(x_pos_point,y_pos_point);
                                h_pos_XY_neg_v2_pT ->Fill(quark_pT);
                            }
                            if(i_quark_mass == 2 && cos_phin > 0.1)
                            {
                                h2D_pos_XY_pos_v2  ->Fill(x_pos_point,y_pos_point);
                                h_pos_XY_pos_v2_pT ->Fill(quark_pT);
                            }

                        }
                        //--------------------------------------------------------------------


                    } // end of event loop
                } // end of quark mass loop
                printf("End of loop \n");
                //--------------------------------------------------------------------

                if(flag_loop)
                {
                    outputfile ->cd();
                    for(Int_t i_mass = 0; i_mass < 8; i_mass++)
                    {
                        cout << "Add profiles to list for i_mass: " << i_mass << endl;
                        tp_v2_vs_pT_mesons[i_mass] ->SetName(Form("v2_vs_pT_BW_id%d_T%d_rho0%d_rhoa%d_Rx%d_fb%d",i_mass,i_Temp,i_rho_0,i_rho_a,i_R_x,i_fboost));
                        list_profiles ->Add((TProfile*)tp_v2_vs_pT_mesons[i_mass]->Clone());
                        //tp_v2_vs_pT_mesons[i_mass] ->Write();
                        tp_v2_vs_pT_mesons[i_mass] ->Reset();

                        h_dN_dpT_mesons[i_mass] ->SetName(Form("h_dN_dpT_vs_pT_BW_id%d_T%d_rho0%d_rhoa%d_Rx%d_fb%d",i_mass,i_Temp,i_rho_0,i_rho_a,i_R_x,i_fboost));
                        list_profiles ->Add((TH1F*)h_dN_dpT_mesons[i_mass]->Clone());
                        //h_dN_dpT_mesons[i_mass] ->Write();
                        h_dN_dpT_mesons[i_mass] ->Reset();
                    }
                }

            }
        }
    }

    if(flag_loop)
    {
        cout << "Write list" << endl;
        outputfile ->cd();
        list_profiles ->Write(Form("list_BW_Rx%d_fb%d",i_R_x,i_fboost),TObject::kSingleKey);
        outputfile ->Close();
        return 1;
    }


    //------------------------------------------------------------------
    // Calculate spectra and v2 analytially
    TGraph* tg_spec = new TGraph();
    TGraph* tg_v2   = new TGraph();
    Double_t inv_yield_BW, v2_BW;
    Double_t integral = 0.0;
    Double_t bin_width = 0.04;
    for(Int_t i_pT = 0; i_pT < 200; i_pT++)
    {
        Double_t pt_BW = i_pT*bin_width;
        blastwave_yield_and_v2(pt_BW, arr_quark_mass_meson[4], Temp_set, rho_0_set, rho_a_set, R_x_set, inv_yield_BW, v2_BW);
        tg_spec ->SetPoint(i_pT,pt_BW,inv_yield_BW*pt_BW);
        tg_v2   ->SetPoint(i_pT,pt_BW,v2_BW);
        integral += bin_width*inv_yield_BW*pt_BW;
    }

    cout << "integral: " << integral << endl;

    for(Int_t i_pT = 0; i_pT < 200; i_pT++)
    {
        Double_t pt_BW;
        tg_spec ->GetPoint(i_pT,pt_BW,inv_yield_BW);
        tg_spec ->SetPoint(i_pT,pt_BW,inv_yield_BW/integral);
    }
    //------------------------------------------------------------------



    h2D_rapidity_vs_eta ->GetXaxis()->SetTitle("#eta");
    h2D_rapidity_vs_eta ->GetYaxis()->SetTitle("y");
    TCanvas* can_rapidity_vs_eta    = Draw_2D_histo_and_canvas(h2D_rapidity_vs_eta,"can_rapidity_vs_eta",1200,700,0,100,"colz");
    PlotLine(-1.0,-1.0,-2.0,2.0,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
    PlotLine(1.0,1.0,-2.0,2.0,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
    PlotLine(-2.0,2.0,0.5,0.5,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
    PlotLine(-2.0,2.0,-0.5,-0.5,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)


    TCanvas* can_single = Draw_1D_histo_and_canvas(tp_v2_vs_pT_mesons[1],"can_single",750,550,0,0.0,"h");

    Int_t plot_particle_dN_dpT = 0;
    h_dN_dpT_mesons[plot_particle_dN_dpT] ->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_dN_dpT_mesons[plot_particle_dN_dpT] ->GetYaxis()->SetTitle("dN/dp_{T} (GeV/c)^{-1}");
    TCanvas* can_dN_dpT_mesons = Draw_1D_histo_and_canvas((TH1D*)h_dN_dpT_mesons[plot_particle_dN_dpT],"can_dN_dpT_mesons",750,550,0,0.0,"h");
    can_dN_dpT_mesons ->cd()->SetLogy(1);
    can_dN_dpT_mesons ->cd();
    h_dN_dpT_jets[plot_particle_dN_dpT] ->SetLineColor(kRed);
    h_dN_dpT_jets[plot_particle_dN_dpT] ->DrawCopy("same h");

    //--------------------------------------------------------------------
    Int_t plot_centrality = 4;
    TH1F* h_dummy = new TH1F("h_dummy","h_dummy",200,0,15);
    h_dummy ->GetXaxis()->SetRangeUser(0.0,15.0);
    h_dummy ->GetYaxis()->SetRangeUser(-0.15,0.45);
    h_dummy ->SetLineColor(10);
    h_dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_dummy->GetYaxis()->SetTitle("v_{2}");
    TCanvas* can_v2_vs_pT = Draw_1D_histo_and_canvas((TH1D*)h_dummy,"can_v2_vs_pT",850,650,0,0.0,"h");

    PlotLine(0.0,15.0,0.0,0.0,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)


    /*
    vec_graphs[plot_centrality] ->SetMarkerColor(arr_color_mass[0]);
    vec_graphs[plot_centrality] ->Draw("same P"); // pions
    //vec_graphs[plot_centrality+7] ->Draw("same P"); // charged kaons
    vec_graphs[plot_centrality+14] ->SetMarkerColor(arr_color_mass[1]);
    vec_graphs[plot_centrality+14] ->Draw("same P"); // K0s
    vec_graphs[plot_centrality+28] ->SetMarkerColor(arr_color_mass[2]);
    vec_graphs[plot_centrality+28] ->Draw("same P"); // protons
    //vec_graphs[plot_centrality+40] ->Draw("same P"); // Lambda
    //vec_graphs[37] ->Draw("same P"); // phi
    //vec_graphs[plot_centrality+47] ->Draw("same P"); // Xi
    //vec_graphs[57] ->Draw("same P"); // Omega
    */
    vec_tge_v2_vs_pT_560_pid[0] ->SetMarkerColor(arr_color_mass[0]);
    //vec_tge_v2_vs_pT_560_pid[0] ->Draw("same P"); // pions
    vec_tge_v2_vs_pT_560_pid[1] ->SetMarkerColor(arr_color_mass[1]);
    //vec_tge_v2_vs_pT_560_pid[1] ->Draw("same P"); // K0s
    vec_tge_v2_vs_pT_560_pid[2] ->SetMarkerColor(arr_color_mass[2]);
    //vec_tge_v2_vs_pT_560_pid[2] ->Draw("same P"); // protons



    for(Int_t i_mass = 0; i_mass < 5; i_mass++)
    {
        //tp_v2_vs_pT_mesons[i_mass] ->Scale(1.6);
        tp_v2_vs_pT_mesons[i_mass] ->SetMarkerColor(arr_color_mass[i_mass]);
        tp_v2_vs_pT_mesons[i_mass] ->SetMarkerStyle(30); // 24
        tp_v2_vs_pT_mesons[i_mass] ->SetLineColor(arr_color_mass[i_mass]);
        tp_v2_vs_pT_mesons[i_mass] ->SetLineStyle(1);
        tp_v2_vs_pT_mesons[i_mass] ->SetLineWidth(2);
        tp_v2_vs_pT_mesons[i_mass] ->DrawCopy("same P");

        Draw_hist_line((TH1D*)tp_v2_vs_pT_mesons[i_mass],0.0,15.0,-0.07,0.7,
                       arr_color_mass[i_mass],4,1,1.0,"ogl");
    }

    tg_JPsi_v2_vs_pT    ->SetMarkerColor(arr_color_mass[3]);
    //tg_JPsi_v2_vs_pT    ->Draw("same P");
    tg_Upsilon_v2_vs_pT ->SetMarkerColor(arr_color_mass[4]);
    //tg_Upsilon_v2_vs_pT ->Draw("same P");

    tg_v2 ->SetLineWidth(4);
    tg_v2 ->SetLineStyle(9);
    tg_v2 ->SetLineColor(kCyan+1);
    tg_v2 ->Draw("same");


    TLegend* leg_v2_vs_pT_A = new TLegend(0.67,0.65,0.77,0.83); // x1,y1,x2,y2
    leg_v2_vs_pT_A->SetBorderSize(0);
    leg_v2_vs_pT_A->SetFillColor(0);
    leg_v2_vs_pT_A->SetTextSize(0.045);
    leg_v2_vs_pT_A->AddEntry((TH1F*)vec_tge_v2_vs_pT_560_pid[0]->Clone(),"#pi","p");
    leg_v2_vs_pT_A->AddEntry((TH1F*)vec_tge_v2_vs_pT_560_pid[1]->Clone(),"K_{s}^{0}","p");
    leg_v2_vs_pT_A->AddEntry((TH1F*)vec_tge_v2_vs_pT_560_pid[2]->Clone(),"p","p");
    leg_v2_vs_pT_A->Draw();

    TLegend* leg_v2_vs_pT_B = new TLegend(0.8,0.71,0.9,0.83); // x1,y1,x2,y2
    leg_v2_vs_pT_B->SetBorderSize(0);
    leg_v2_vs_pT_B->SetFillColor(0);
    leg_v2_vs_pT_B->SetTextSize(0.045);
    leg_v2_vs_pT_B->AddEntry((TGraphErrors*)tg_Upsilon_v2_vs_pT->Clone(),"#Upsilon","p");
    leg_v2_vs_pT_B->AddEntry((TGraphErrors*)tg_JPsi_v2_vs_pT->Clone(),"J/#Psi","p");
    leg_v2_vs_pT_B->Draw();


    plotTopLegend((char*)"|y|<0.5",0.66,0.83,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    plotTopLegend((char*)"2.5<y<4",0.79,0.83,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    plotTopLegend((char*)"5-60%",0.25,0.83,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1



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

#if 0
    if(tp_v2_vs_pT_opt_best)
    {
        printf("Plot tp_v2_vs_pT_opt_best \n");
        tp_v2_vs_pT_opt_best  ->SetMarkerColor(kGreen+1);
        tp_v2_vs_pT_opt_best  ->SetMarkerStyle(29);
        tp_v2_vs_pT_opt_best  ->SetMarkerSize(1.2);
        tp_v2_vs_pT_opt_best  ->Draw("same P");
    }
#endif
    //--------------------------------------------------------------------


    //--------------------------------------------------------------------
    plot_spectra();
    //--------------------------------------------------------------------


#if 0
    //--------------------------------------------------------------------
    TCanvas* can_v2_vs_pT_compare = Draw_1D_graph_and_canvas(vec_tg_v2_vs_pT_Mathematica[0],"can_v2_vs_pT_compare",800,800,-0.05,0.42,"",20,0.9,kBlack);

    for(Int_t i_mass = 0; i_mass < 3; i_mass++)
    {
        vec_tg_v2_vs_pT_Mathematica[i_mass] ->SetMarkerColor(arr_color_mass[i_mass]);
        vec_tg_v2_vs_pT_Mathematica[i_mass] ->SetMarkerStyle(28); // 24
        vec_tg_v2_vs_pT_Mathematica[i_mass] ->SetLineColor(arr_color_mass[i_mass]);
        vec_tg_v2_vs_pT_Mathematica[i_mass] ->SetLineStyle(1);
        vec_tg_v2_vs_pT_Mathematica[i_mass] ->SetLineWidth(2);
        vec_tg_v2_vs_pT_Mathematica[i_mass] ->Draw("same P");
    }

    for(Int_t i_mass = 0; i_mass < 3; i_mass++)
    {
        //tp_v2_vs_pT_mesons[i_mass] ->Scale(1.6);
        tp_v2_vs_pT_mesons[i_mass] ->SetMarkerColor(arr_color_mass[i_mass]);
        tp_v2_vs_pT_mesons[i_mass] ->SetMarkerStyle(30); // 24
        tp_v2_vs_pT_mesons[i_mass] ->SetLineColor(arr_color_mass[i_mass]);
        tp_v2_vs_pT_mesons[i_mass] ->SetLineStyle(1);
        tp_v2_vs_pT_mesons[i_mass] ->SetLineWidth(2);
        tp_v2_vs_pT_mesons[i_mass] ->Draw("same P");
    }
    //--------------------------------------------------------------------
#endif



#if 0
    //--------------------------------------------------------------------
    TCanvas* can_v2_vs_pt = new TCanvas("can_v2_vs_pt","can_v2_vs_pt",700,10,500,500);
    can_v2_vs_pt ->cd();
    for(Int_t i_mass = 0; i_mass < N_masses; i_mass++)
    {
        tp_v2_vs_pT_mesons[i_mass] ->SetMarkerColor(arr_color_mass[i_mass]);
        tp_v2_vs_pT_mesons[i_mass] ->SetMarkerStyle(20);
        tp_v2_vs_pT_mesons[i_mass] ->Draw("same P");
    }
    //--------------------------------------------------------------------
#endif


#if 0
    //--------------------------------------------------------------------
    TCanvas* can_dN_dphi_vs_phi = new TCanvas("can_dN_dphi_vs_phi","can_dN_dphi_vs_phi",700,10,500,500);
    can_dN_dphi_vs_phi ->cd();
    tp_dN_dphi_vs_phi ->Scale(1.0/tp_dN_dphi_vs_phi ->GetEntries());
    tp_dN_dphi_vs_phi ->DrawCopy();

    for(Int_t i_par = 0; i_par < 5; i_par++)
    {
        f_FlowFitFunc->ReleaseParameter(i_par);
        f_FlowFitFunc->SetParError(i_par,0.0);
        f_FlowFitFunc->SetParameter(i_par,0.0);
    }

    f_FlowFitFunc ->SetParameter(0,1.0);
    f_FlowFitFunc ->FixParameter(1,0.0);
    f_FlowFitFunc ->SetParameter(2,-0.7);
    f_FlowFitFunc ->FixParameter(3,0.0);
    f_FlowFitFunc ->SetParameter(4,0.0);

    f_FlowFitFunc ->SetRange(-TMath::Pi(),TMath::Pi());
    tp_dN_dphi_vs_phi ->Fit("f_FlowFitFunc","WMN","",-TMath::Pi(),TMath::Pi());
    Double_t v2_fit  = f_FlowFitFunc->GetParameter(1);

    f_FlowFitFunc ->SetLineColor(kGreen);
    f_FlowFitFunc ->SetLineStyle(2);
    f_FlowFitFunc ->SetLineWidth(2);
    f_FlowFitFunc ->DrawCopy("same");
    //--------------------------------------------------------------------



    //--------------------------------------------------------------------
    TCanvas* can_pos_XY_neg_v2    = Draw_2D_histo_and_canvas(h2D_pos_XY_neg_v2,"can_pos_XY_neg_v2",1600,1000,0,100,"colz");
    TCanvas* can_pos_XY_neg_v2_pT = Draw_1D_histo_and_canvas(h_pos_XY_neg_v2_pT,"can_pos_XY_neg_v2_pT",750,550,0,0.0,"h");
    //--------------------------------------------------------------------


    //--------------------------------------------------------------------
    TCanvas* can_pos_XY_pos_v2    = Draw_2D_histo_and_canvas(h2D_pos_XY_pos_v2,"can_pos_XY_pos_v2",1600,1000,0,100,"colz");
    TCanvas* can_pos_XY_pos_v2_pT = Draw_1D_histo_and_canvas(h_pos_XY_pos_v2_pT,"can_pos_XY_pos_v2_pT",750,550,0,0.0,"h");
    //--------------------------------------------------------------------
#endif


#if 0
    //--------------------------------------------------------------------
    TCanvas*  can_cos_theta_vs_pT   = Draw_2D_histo_and_canvas(h2D_cos_theta_vs_pT,"can_cos_theta_vs_pT",1600,1000,0,100,"colz");
    //--------------------------------------------------------------------
#endif

    outputfile ->cd();
    for(Int_t i_mass = 0; i_mass < 5; i_mass++)
    {
        tp_v2_vs_pT_mesons[i_mass] ->SetName(Form("v2_vs_pT_BW_id%d_T%4.2f",i_mass,Temp_set));
        tp_v2_vs_pT_mesons[i_mass] ->Write();
    }
    outputfile ->Close();

    return 0;

}