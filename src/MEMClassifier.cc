#include "TTH/CommonClassifier/interface/MEMClassifier.h"

MEMResult MEMClassifier::GetOutput(
    const std::vector<TLorentzVector>& selectedLeptonP4,
    const std::vector<TLorentzVector>& selectedJetP4,
    const std::vector<double>& selectedJetCSV,
    const std::vector<TLorentzVector>& looseSelectedJetP4,
    const std::vector<double>& looseSelectedJetCSV,
    TLorentzVector& metP4
) {
    //integrand->next_event();
    if (selectedLeptonP4.size() != 1) {
        throw std::runtime_error("Expected a single-lepton event");
    }
    
    //store MEM objects
    std::vector<MEM::Object*> objs;

    for (unsigned int ij=0; ij<selectedJetP4.size(); ij++) {
        TLorentzVector p4 = selectedJetP4.at(ij);
        assert(p4.Pt() > 0);
        MEM::Object* jet = make_jet(
            p4.Pt(), p4.Eta(), p4.Phi(), p4.M(), selectedJetCSV.at(ij)
        );
        objs.push_back(jet);
        integrand->push_back_object(jet);
    }

    for (auto& lep_p4 : selectedLeptonP4) {
        assert(lep_p4.Pt() > 0);
        MEM::Object* lep = make_lepton(lep_p4.Pt(), lep_p4.Eta(), lep_p4.Phi(), lep_p4.M(), +1);
        objs.push_back(lep);
        integrand->push_back_object(lep);
    }

    assert(metP4.Pt() > 0);
    MEM::Object met(metP4, MEM::ObjectType::MET );
    integrand->push_back_object(&met);
    integrand->set_cfg(cfg);
    cout << "cfg.transfer_function_method=" << cfg.transfer_function_method << endl;
    integrand->set_permutation_strategy
    ({MEM::Permutations::BTagged,
      MEM::Permutations::QUntagged,
      MEM::Permutations::QQbarBBbarSymmetry
    });

    std::cout << "evaluating signal" << std::endl;
    MEM::MEMOutput res_sig = integrand->run( MEM::FinalState::LH, MEM::Hypothesis::TTH, {});
    std::cout << "evaluating bkg" << std::endl;
    MEM::MEMOutput res_bkg = integrand->run( MEM::FinalState::LH, MEM::Hypothesis::TTBB, {});

    MEMResult res;
    res.p = res_sig.p / (res_sig.p + 0.02*res_bkg.p);
    return res;
}
  
MEM::Object* MEMClassifier::make_jet(double pt, double eta, double phi, double mass, double csv) const {
    TLorentzVector lv;
    lv.SetPtEtaPhiM(pt, eta, phi, mass);
    MEM::Object* obj = new MEM::Object(lv, MEM::ObjectType::Jet);
    obj->addObs( MEM::Observable::BTAG, csv); // 0 - jet is assumed to be from a light quark, 1 - a b quark
    obj->addObs( MEM::Observable::PDGID, 0);  // currently not used
    // attach the transfer functions corresponding to the jet
    obj->addTransferFunction(MEM::TFType::bReco, getTransferFunction("b", lv.Eta()));
    obj->addTransferFunction(MEM::TFType::qReco, getTransferFunction("l", lv.Eta()));
    return obj;
}

MEM::Object* MEMClassifier::make_lepton(double pt, double eta, double phi, double mass, double charge) const {
    TLorentzVector lv;
    lv.SetPtEtaPhiM(pt, eta, phi, mass);
    MEM::Object* obj = new MEM::Object(lv, MEM::ObjectType::Lepton );
    obj->addObs( MEM::Observable::CHARGE, charge);
    return obj;
}

// Returns the transfer function corresponding to a jet flavour and eta
TF1* MEMClassifier::getTransferFunction(const char* flavour, double eta) const {
    int etabin = 0;
    if (std::abs(eta) > 1.0) {
        etabin = 1;
    }
    stringstream ss;
    ss << "tf_" << flavour << "_etabin" << etabin;
    const char* fname = ss.str().c_str();
    TF1* tf = (TF1*)(transfers->Get(fname));
    if (tf == 0) {
        cerr << "could not get transfer function " << fname << endl;
        cerr << flush;
        throw std::exception();
    }
    return tf;
}

MEMClassifier::MEMClassifier() : cfg(MEM::MEMConfig()) {
    cout << "Calling MEMClassifier()" << endl;

    transfers = TFile::Open("root/transfer.root");
    cfg.defaultCfg();
    cfg.transfer_function_method = MEM::TFMethod::External;
    //Transfer functions for jet reconstruction efficiency
    cfg.set_tf_global(MEM::TFType::bLost, 0, getTransferFunction("beff", 0.0));
    cfg.set_tf_global(MEM::TFType::bLost, 1, getTransferFunction("beff", 2.0));
    cfg.set_tf_global(MEM::TFType::qLost, 0, getTransferFunction("leff", 0.0));
    cfg.set_tf_global(MEM::TFType::qLost, 1, getTransferFunction("leff", 2.0));

    integrand = new MEM::Integrand( 
        MEM::DebugVerbosity::output
        //|MEM::DebugVerbosity::init
        //|MEM::DebugVerbosity::input
        //|MEM::DebugVerbosity::init_more
        //|MEM::DebugVerbosity::integration,
        ,cfg
    );
    integrand->set_cfg(cfg);

}

MEMClassifier::~MEMClassifier() {
    delete integrand;
    delete transfers;
}
