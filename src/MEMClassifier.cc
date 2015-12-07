#include "TTH/CommonClassifier/interface/MEMClassifier.h"

MEMClassifier::GetOutput(
    const std::vector<TLorentzVector>& selectedLeptonP4,
    const std::vector<TLorentzVector>& selectedJetP4,
    const std::vector<double>& selectedJetCSV,
    const std::vector<TLorentzVector>& looseSelectedJetP4,
    const std::vector<double>& looseSelectedJetCSV,
    TLorentzVector& metP4
) {
    if (selectedLeptonP4.size() != 1) {
        throw std::exception("Expected a single-lepton event");
    }

    for (auto& lep_p4 : selectedLeptonP4) {
    }
}
  
Object MEMClassifier::make_jet(double pt, double eta, double phi, double mass, double csv, TFile* tffile) const {
    TLorentzVector lv_j1;
    lv_j1.SetPtEtaPhiM(pt, eta, phi, mass);
    Object j1( lv_j1, ObjectType::Jet );
    j1.addObs( Observable::BTAG, csv); // 0 - jet is assumed to be from a light quark, 1 - a b quark
    j1.addObs( Observable::PDGID, 0);  // currently not used
    // attach the transfer functions corresponding to the jet
    j1.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j1.Eta()));
    j1.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j1.Eta()));
    return j1;
}

Object MEMClassifier::make_lepton(double pt, double eta, double phi, double mass, double charge, TFile* tffile) {
    TLorentzVector lv;
    lv.SetPtEtaPhiM(pt, eta, phi, mass);
    Object l(lv, ObjectType::Lepton );
    l.addObs( Observable::CHARGE, charge);
    return l;
}

// Returns the transfer function corresponding to a jet flavour and eta
TF1* MEMClassifier::getTransferFunction(TFile* tffile, const char* flavour, double eta) const {
    int etabin = 0;
    if (std::abs(eta) > 1.0) {
        etabin = 1;
    }
    stringstream ss;
    ss << "tf_" << flavour << "_etabin" << etabin;
    const char* fname = ss.str().c_str();
    TF1* tf = (TF1*)(tffile->Get(fname));
    if (tf == 0) {
        cerr << "could not get transfer function " << fname << endl;
        cerr << flush;
        throw exception();
    }
    return tf;
}
