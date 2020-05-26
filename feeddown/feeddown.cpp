#include "feeddown.h"

// ClassImp(feeddown)

//
// feeddown class: see feeddown.h for a general description
//

void feeddown::set_model_parameters(Double_t p[], const Int_t npar) {
    if (npar <= 20) {
        // copy model parameters
        for (Int_t i = 0; i < npar; ++i)
            mp[i] = p[i];
    } else {
        std::cout << "ERROR: too many model parameters" << std::endl;
    }
}

// set model for calculating the yield (dN/dpT)
void feeddown::set_dndpt_model(Double_t model(Double_t *x, Double_t *par)) {
    dndpt_model = model;
    current_model = dndpt;
}

// set model for simultaneously calculating yield (yield dN/dpT and v2)
void feeddown::set_dndpt_and_v2_model(void model(Double_t *x, Double_t *par, Double_t &inv_yield, Double_t &v2)) {
    dndpt_and_v2_model = model;
    current_model = dndpt_and_v2;
}

// Fill dn/dpT histogram for primary particles for the given model
void feeddown::set_dndpt_and_v2_primary_from_model(TH1D *h_dndpt, TH1D *h_v2) {

    if (current_model == dndpt) {
        for (Int_t i = 1; i <= h_dndpt->GetXaxis()->GetNbins(); ++i) {
            Double_t pt = h_dndpt->GetBinCenter(i);
            Double_t yield = (*dndpt_model)(&pt, mp);
            h_dndpt->SetBinContent(i, yield);
            h_v2->SetBinContent(i, 0);
        }
    } 
    else if (current_model == dndpt_and_v2) {
        Double_t dndpt = 0;
        Double_t v2 = 0;
        for (Int_t i = 1; i <= h_dndpt->GetXaxis()->GetNbins(); ++i) {
            Double_t pt = h_dndpt->GetBinCenter(i);
            (*dndpt_and_v2_model)(&pt, mp, dndpt, v2);
            h_dndpt->SetBinContent(i, dndpt);
            h_v2->SetBinContent(i, v2);
        }

    } else {
        std::cout << "ERROR: model not defined" << std::endl;
    }
}

// read 2d decay histograms from root file
void feeddown::read_decay_histograms(TString filename) {

    TFile f(filename);
    TIter next(f.GetListOfKeys());
    TKey *key;

    while ((key = (TKey *)next())) {

        std::string class_name = key->GetClassName();
        std::string name = key->GetName();

        // check for 2d decay histograms
        const std::regex pieces_regex_2d("decay_hist_dndpt_secondary_idprim([-]*[1-9]+)_idsec([-]*[1-9]+)");
        std::smatch pieces_match_2d;
        if (std::regex_match(name, pieces_match_2d, pieces_regex_2d)) {
            Int_t idprim = std::stoi(pieces_match_2d[1].str());
            Int_t idsec = std::stoi(pieces_match_2d[2].str());
            decay_hist_dndpt_secondary[idprim][idsec] = (TH2D *)key->ReadObj();
            decay_hist_dndpt_secondary[idprim][idsec]->SetDirectory(0); // decouple histo from the open file director
        }

        // check for 2d cos2phi histograms
        const std::regex pieces_regex_v2("decay_hist_cos2phi_secondary_idprim([-]*[1-9]+)_idsec([-]*[1-9]+)");
        std::smatch pieces_match_v2;
        if (std::regex_match(name, pieces_match_v2, pieces_regex_v2)) {
            Int_t idprim = std::stoi(pieces_match_v2[1].str());
            Int_t idsec = std::stoi(pieces_match_v2[2].str());
            decay_hist_cos2phi_secondary[idprim][idsec] = (TH2D *)key->ReadObj();
            decay_hist_cos2phi_secondary[idprim][idsec]->SetDirectory(0); // decouple histo from the open file director
        }

        // get 1d normalization histograms
        const std::regex pieces_regex_1d("decay_hist_dndpt_primary_id([-]*[1-9]+)");
        std::smatch pieces_match_1d;
        if (std::regex_match(name, pieces_match_1d, pieces_regex_1d)) {
            Int_t id = std::stoi(pieces_match_1d[1].str());
            decay_hist_dndpt_primary[id] = (TH1D *)key->ReadObj();
            decay_hist_dndpt_primary[id]->SetDirectory(0);
        }
    }

    f.Close();

    // copy histograms to get the correct number of bins and histogram limits
    // pick a particle id that should be always present (pi0 = 111) 
    h_dndpt_primary = *decay_hist_dndpt_primary[111];
    h_dndpt_secondary = *decay_hist_dndpt_primary[111];
    h_dndpt_total = *decay_hist_dndpt_primary[111];

    h_v2_primary = *decay_hist_dndpt_primary[111];
    h_v2_secondary = *decay_hist_dndpt_primary[111];
    h_v2_total = *decay_hist_dndpt_primary[111];
}

// calculate feeddown contribution for particle species id
// depending on the model type (dndpt or dndpt_and_v2) either the feeddown
// only for dndpt or dndpt and v2 is calculated
void feeddown::calc_feeddown_hist(const Int_t id) {

    TString particle_name = pythia->Pythia8()->particleData.name(id);
    TString title_prim = particle_name + " primary";
    TString title_sec = particle_name + " secondary";

    // reset hisograms for dndpt and v2
    h_dndpt_primary.Reset();
    h_dndpt_secondary.Reset();
    h_dndpt_total.Reset();
    h_v2_primary.Reset();
    h_v2_secondary.Reset();
    h_v2_total.Reset();

    //  set histograms names
    h_dndpt_primary.SetTitle(title_prim.Data());
    h_dndpt_secondary.SetTitle(title_sec.Data());
    h_v2_primary.SetTitle(title_prim.Data());
    h_v2_secondary.SetTitle(title_sec.Data());

    // mass and spin type of the considered primary particle
    Double_t m = pythia->Pythia8()->particleData.m0(id);
    Double_t spinType = pythia->Pythia8()->particleData.spinType(id); // spinType = 2 J + 1
    mp[0] = m;
    mp[1] = spinType;
    set_dndpt_and_v2_primary_from_model(&h_dndpt_primary, &h_v2_primary);

    const Int_t n_bins = h_dndpt_secondary.GetNbinsX();
    // Double_t n_sum_mother = 0; // mother particle yields summed for all species and pT bins
    Double_t dndpt_sum_sec[n_bins + 1]; // secondary yield for differnt pT,sec bins
    Double_t v2_sum_sec[n_bins + 1];    // secondary v2 sum for differnt pT,sec bins

    calc_feeddown_in_pt_sec_index_range(id, 1, n_bins, dndpt_sum_sec, v2_sum_sec);

    // set class member variables (sum of pt-integrated yields for all mother particles)
    // n_secondary = n_sum_mother;
    // n_total = n_primary + n_secondary;

    // store dndpt (and v2) results into the corresponding histograms
    for (Int_t j = 1; j <= n_bins; ++j) {

        // dn/dpt secondary
        h_dndpt_secondary.SetBinContent(j, dndpt_sum_sec[j]);

        // dn/dpt total
        Double_t dndpt_primary = h_dndpt_primary.GetBinContent(j);
        h_dndpt_total.SetBinContent(j, dndpt_primary + dndpt_sum_sec[j]);
        // h_dndpt_total->SetBinContent(j, 2. * dndpt_primary);

        // v2 secondary
        if (current_model == dndpt_and_v2) {
            Double_t v2_secondary = 0;
            if (dndpt_sum_sec[j] != 0)
                v2_secondary = v2_sum_sec[j] / dndpt_sum_sec[j];
            h_v2_secondary.SetBinContent(j, v2_secondary);

            // v2 total
            Double_t v2_primary = h_v2_primary.GetBinContent(j);
            Double_t v2_total =
                (dndpt_primary * v2_primary + dndpt_sum_sec[j] * v2_secondary) / (dndpt_primary + dndpt_sum_sec[j]);
            h_v2_total.SetBinContent(j, v2_total);
        }
    }
}

// feeddown calculation for fixed pT
// interface for dn/dpt variables only
void feeddown::calc_feeddown(const Int_t id, const Double_t &pt, Double_t &dndpt_primary, Double_t &dndpt_secondary,
                             Double_t &dndpt_total) {

    Double_t v2_primary_dummy = 0;
    Double_t v2_secondary_dummy = 0;
    Double_t v2_total_dummy = 0;

    calc_feeddown(id, pt, dndpt_primary, dndpt_secondary, dndpt_total, v2_primary_dummy, v2_secondary_dummy,
                  v2_total_dummy);
};

// feeddown calculation for fixed pT
// interface for dn/dpt and v2 variables
void feeddown::calc_feeddown(const Int_t id, const Double_t &pt, Double_t &dndpt_primary, Double_t &dndpt_secondary,
                             Double_t &dndpt_total, Double_t &v2_primary, Double_t &v2_secondary, Double_t &v2_total) {

    // mass and spin type of the considered primary particle
    Double_t m = pythia->Pythia8()->particleData.m0(id);
    Double_t spinType = pythia->Pythia8()->particleData.spinType(id); // spinType = 2 J + 1
    mp[0] = m;
    mp[1] = spinType;

    dndpt_primary = 0;
    v2_primary = 0;

    Double_t pttmp = pt;
    if (current_model == dndpt) {
        dndpt_primary = (*dndpt_model)(&pttmp, mp);
    } else if (current_model == dndpt_and_v2) {
        (*dndpt_and_v2_model)(&pttmp, mp, dndpt_primary, v2_primary);
    } else {
        std::cout << "ERROR: model not defined" << std::endl;
    }

    // Double_t n_sum_mother = 0; // mother particle yields summed for all species and pT bins
    const Int_t n_bins = h_dndpt_secondary.GetNbinsX();
    Double_t dndpt_sum_sec[n_bins + 1]; // secondary yield for differnt pT,sec bins
    Double_t v2_sum_sec[n_bins + 1];    // secondary v2 sum for differnt pT,sec bins

    const Int_t j_bin = h_dndpt_primary.FindBin(pt);
    calc_feeddown_in_pt_sec_index_range(id, j_bin, j_bin, dndpt_sum_sec, v2_sum_sec);

    // dn/dpT
    dndpt_secondary = dndpt_sum_sec[j_bin];
    dndpt_total = dndpt_primary + dndpt_secondary;

    if (current_model == dndpt_and_v2) {
        v2_secondary = 0;
        if (dndpt_sum_sec[j_bin] != 0)
            v2_secondary = v2_sum_sec[j_bin] / dndpt_sum_sec[j_bin];
        v2_total = (dndpt_primary * v2_primary + dndpt_secondary * v2_secondary) / (dndpt_primary + dndpt_secondary);
    }
};

void feeddown::calc_feeddown_in_pt_sec_index_range(const Int_t id, const Int_t j_min, Int_t j_max,
                                                   Double_t dndpt_sum_sec[], Double_t v2_sum_sec[]) {

    // histogram to store pt and v2 for the mother particle
    TH1D *h_dndpt_mother = (TH1D *)decay_hist_dndpt_primary[111]->Clone();
    TH1D *h_v2_mother = (TH1D *)decay_hist_dndpt_primary[111]->Clone();

    const Int_t n_bins = h_dndpt_mother->GetNbinsX();

    // set counters to zero (index 0 is not used as we are dealing with root histograms)
    for (Int_t j = j_min; j <= j_max; ++j) {
        dndpt_sum_sec[j] = 0;
        v2_sum_sec[j] = 0;
    }

    // loop over available mother id's
    for (auto const &mothers : decay_hist_dndpt_secondary) {

        // id, mass, and spin of the mother particles
        Int_t id_mother = mothers.first;
        Double_t m = pythia->Pythia8()->particleData.m0(id_mother);
        Double_t spinType = pythia->Pythia8()->particleData.spinType(id_mother); // spinType = 2 J + 1
  
        // reject particle if its mass is above the set mass limts
  		if (m > maximum_mother_mass) continue; 

        // if set of mother particle is not empty check whether mother particle is in the set
    	if (mother_particles.size() != 0 && mother_particles.count(id_mother) == 0) continue;
     
        mp[0] = m;
        mp[1] = spinType;

        // calculate dndpt and v2 for this mother particle
        h_dndpt_mother->Reset();
        h_v2_mother->Reset();
        set_dndpt_and_v2_primary_from_model(h_dndpt_mother, h_v2_mother);
        // n_sum_mother += n_mother;

        // add feeddown from species id_mother
        if (id_mother != id && decay_hist_dndpt_secondary[id_mother].count(id) > 0) {

            // i = index of mother pT bin
            // j = index of daughter pT bin

            // loop over mother pT
            for (Int_t i = 1; i <= n_bins; ++i) {

                Double_t dndpt_mother = h_dndpt_mother->GetBinContent(i);
                Double_t v2_mother = 0;
                if (current_model == dndpt_and_v2)
                    v2_mother = h_v2_mother->GetBinContent(i);

                // loop over daughter
                // for heavy particles the daughter can have higher pT than mother, therefore we start a index 0
                for (Int_t j = j_min; j <= j_max; ++j) {

                    // dndpt
                    Double_t dndpt_sec = decay_hist_dndpt_secondary[id_mother][id]->GetBinContent(i, j);
                    if (dndpt_sec == 0)
                        continue;

                    Double_t dndpt_sec_norm = dndpt_sec / decay_hist_dndpt_primary[id_mother]->GetBinContent(i);
                    dndpt_sum_sec[j] += dndpt_mother * dndpt_sec_norm;

                    // v2 calculation in case a v2 model is defined
                    if (current_model == dndpt_and_v2) {
                        Double_t cos2phi_sum = decay_hist_cos2phi_secondary[id_mother][id]->GetBinContent(i, j);
                        if (cos2phi_sum == 0)
                        	continue;

                        Double_t v2delta = 0;
                        if (dndpt_sec != 0)
                            v2delta = cos2phi_sum / dndpt_sec;

                        v2_sum_sec[j] += dndpt_mother * dndpt_sec_norm * v2_mother * v2delta;
                    }

                } // loop over daughter pT
            }
        } // if

        // debug
        // cout << id_mother << " " << h_dndpt_secondary->GetBinContent(10) << endl;
    } // loop over mother particles

    delete h_dndpt_mother;
    delete h_v2_mother;
}
