// Class for the calculation of the feeddown contribution to transverse momentum spectra and
// elliptic flow v2(pt). The class first reads in two-dimensional decay histograms. A 2d decay histogram
// describes the pT spectra of hadron species Y resulting from electromagnetic and strong decays of
// the mother species X. The total feeddown contribution to species y is the sum over the contribution of all mother
// species.
// 

#pragma clang diagnostic ignored "-Wavailability"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TPythia8.h"
#include <regex>
#include <string>

class feeddown {

    TPythia8 *pythia;

    // dn/dpT and/or dn/dpT + v2 model
    // parameter 0 is assumed to be the particle mass
    // parameter 1 is assumed to be the spin type = 2 S + 1
    Double_t mp[20]; // array with model parameters, at maximum 20 parameters

    // models need to have a certain format (see below)
    Double_t (*dndpt_model)(Double_t *x, Double_t *par);
    void (*dndpt_and_v2_model)(Double_t *x, Double_t *par, Double_t &dndpt, Double_t &v2);

    enum model_type { undefined = 0, dndpt = 1, dndpt_and_v2 = 2 };
    model_type current_model;

    // two-dimensional decay histograms:
    // x axis of the decay histogram: pT of the mother particle X. Here a flat
    // pT distribution is used. y axis of the decay histogram: pT of the daughter particle.
    std::map<Int_t, TH1D *> decay_hist_dndpt_primary;
    std::map<Int_t, std::map<Int_t, TH2D *>> decay_hist_dndpt_secondary;
    std::map<Int_t, std::map<Int_t, TH2D *>> decay_hist_cos2phi_secondary;

    // list of considered mother particles
    // if empty, all mother particles in the decay histogram file are considered 
    std::unordered_set<Int_t> mother_particles;

    // possibilty to set upper mass limit for considered mother particles
    Double_t maximum_mother_mass;

    // for a given model, fill histogram for dndpt (and v2 if v2 model is defined)
    // for the current parameter set (as given by mp)
    void set_dndpt_and_v2_primary_from_model(TH1D *, TH1D *);
    void calc_feeddown_in_pt_sec_index_range(const Int_t, const Int_t, Int_t, Double_t *, Double_t *);

    TH1D h_dndpt_primary;   // dn/dpT_(primary_particles)
    TH1D h_dndpt_secondary; // dn/dpT_(secondary_particles)
    TH1D h_dndpt_total;     // dn/dpT_(total) = dn/dpT_(primary_particles) + dn/dpT_(secondary_particles)

    TH1D h_v2_primary;   // v2sum (primary_particles)
    TH1D h_v2_secondary; // v2sum (secondary_particles)
    TH1D h_v2_total;     //

  public:
    feeddown() : pythia{new TPythia8}, current_model{undefined}, maximum_mother_mass{10000.}{};

    void read_decay_histograms(TString filename);

    // feeddown for the full pt range
    void calc_feeddown_hist(const Int_t);

    // feeddown for a given pt bin (only yield dn/dpT or both)
    void calc_feeddown(const Int_t, const Double_t &pt, Double_t &dndpt_primary, Double_t &dndpt_secondary,
                       Double_t &dndpt_total);
    void calc_feeddown(const Int_t, const Double_t &pt, Double_t &dndpt_primary, Double_t &dndpt_secondary,
                       Double_t &dndpt_total, Double_t &v2_primary, Double_t &v2_secondary, Double_t &v2_total);

    void set_dndpt_model(Double_t model(Double_t *x, Double_t *par));
    void set_dndpt_and_v2_model(void model(Double_t *x, Double_t *par, Double_t &dndpt, Double_t &v2));
    void set_model_parameters(Double_t p[], const Int_t npar);

    // possibility to provide list of considered mother particles
    // if nothing is provided, the full list of mother particles in the decay histogram file is consideres
    void set_mother_particles(std::unordered_set<Int_t> mothers) {mother_particles = mothers;};
    void reset_mother_particles() {mother_particles.clear();};

    void set_maximum_mother_mass(Double_t m_max) {maximum_mother_mass = m_max;};
    void reset_maximum_mother_mass() {maximum_mother_mass = 10000;}; // mass = 10000 GeV --> no particle rejected

    // methods for retrieving the results from histograms
    // histograms pointers
    const TH1D *get_dndpt_primary_hist() { return &h_dndpt_primary; };
    const TH1D *get_dndpt_secondary_hist() { return &h_dndpt_secondary; };
    const TH1D *get_dndpt_total_hist() { return &h_dndpt_total; };
    const TH1D *get_v2_primary_hist() { return &h_v2_primary; };
    const TH1D *get_v2_secondary_hist() { return &h_v2_secondary; };
    const TH1D *get_v2_total_hist() { return &h_v2_total; };

    // values at a specific transverse momentum from calculated histograms
    Double_t get_dndpt_primary_from_hist(const Double_t &pt) {return h_dndpt_primary.Interpolate(pt);};
    Double_t get_dndpt_secondary_from_hist(const Double_t &pt) {return h_dndpt_secondary.Interpolate(pt); };
    Double_t get_dndpt_total_from_hist(const Double_t &pt) {return h_dndpt_total.Interpolate(pt);};
    Double_t get_v2_primary_from_hist(const Double_t &pt) {return h_v2_primary.Interpolate(pt);};
    Double_t get_v2_secondary_from_hist(const Double_t &pt) {return h_v2_secondary.Interpolate(pt); };
    Double_t get_v2_total_from_hist(const Double_t &pt) {return h_v2_total.Interpolate(pt);};

    // ClassDef(feeddown,1);
};