#include "TTH/CommonClassifier/interface/MEMClassifier.h"
#include "TFile.h"
#include <iostream>

using namespace std;
using namespace MEM;

int main(){

  //Load the transfer functions
  TFile* tffile = new TFile("root/transfer.root");

  //create a MEM configuration.
  //this needs to be done once per job, not for every event
  MEMConfig cfg;
  cfg.defaultCfg();
  cfg.transfer_function_method = TFMethod::External;

  //Transfer functions for jet reconstruction efficiency
  cfg.set_tf_global(TFType::bLost, 0, *getTransferFunction(tffile, "beff", 0.0));
  cfg.set_tf_global(TFType::bLost, 1, *getTransferFunction(tffile, "beff", 2.0));
  cfg.set_tf_global(TFType::qLost, 0, *getTransferFunction(tffile, "leff", 0.0));
  cfg.set_tf_global(TFType::qLost, 1, *getTransferFunction(tffile, "leff", 2.0));

  //Create the mem integrator, once per job
  Integrand* integrand = new Integrand( 
    DebugVerbosity::output
    //|DebugVerbosity::init
    //|DebugVerbosity::input
    //|DebugVerbosity::init_more
    //|DebugVerbosity::integration				        
    ,cfg
  );

  //Add some objects to the MEM
  Object j1 = make_jet(242.816604614, -0.107542805374, 1.25506973267, 24.5408706665, 0.0, tffile);
  Object j2 = make_jet(191.423553467, -0.46368226409, 0.750520706177, 30.5682048798, 0.0, tffile);
  Object j3 = make_jet(77.6708831787, -0.709680855274, -2.53739523888, 10.4904966354, 1.0, tffile);
  Object j4 = make_jet(235.892044067, -0.997860729694, -2.10646605492, 27.9887943268, 1.0, tffile);
  Object j5 = make_jet(52.0134391785, -0.617823541164, -1.23360788822, 6.45914268494, 1.0, tffile);
  Object j6 = make_jet(35.6511192322, 0.566395223141, -2.51394343376, 8.94268417358, 1.0, tffile);
  
  Object l1 = make_lepton(52.8751449585, -0.260020583868, -2.55171084404, 0.139569997787, +1.0, tffile);

  //create a MET
  TLorentzVector lv_met;
  lv_met.SetPtEtaPhiM( 92.1731872559,0., -1.08158898354, 0.);
  Object met( lv_met, ObjectType::MET );
  
  //add all objects to the MEM integrator
  integrand->push_back_object( &j1 );
  integrand->push_back_object( &j2 );
  integrand->push_back_object( &j3 );
  integrand->push_back_object( &j4 );
  integrand->push_back_object( &j5 );
  integrand->push_back_object( &j6 );
  integrand->push_back_object( &l1 );
  integrand->push_back_object( &met );
  
  //permute light quarks only in untagged jets, b-quarks in tagged
  //remove symmetric permutations
  integrand->set_permutation_strategy
    (  {Permutations::BTagged
        ,Permutations::QUntagged 
        ,Permutations::QQbarBBbarSymmetry
        } 
    );

  MEMOutput res;			   
  
  //Evaluate fully reconstructed hypothesis
  //LH - single-leptonic decay channel
  //LL - dileptonic decay channel
  cout << "Fully reconstructed interpretation" << endl;
  cout << "evaluating tth hypo" << endl;
  //third variable is variables to integrate over
  //if nothing is specified, assume that all jets (4b + 2 light) were reconstructed
  res = integrand->run( FinalState::LH, Hypothesis::TTH,  {} );
  cout << "p = " << res.p << " +- " << res.p_err << endl;
  double p0 = res.p;

  cout << "evaluating ttbb hypo" << endl;
  res = integrand->run( FinalState::LH, Hypothesis::TTBB, {} );
  cout << "p = " << res.p << " +- " << res.p_err << endl;
  double p1 = res.p;

  //this is the final discriminator value. the normalization constant is to be optimized
  double mem_w = p0 / (p0 + 0.02*p1);
  cout << "mem 222 discriminator " << mem_w << endl;

  ////Evaluate 022 hypothesis. We do not use the information provided by the light quarks.
  //cout << "Integrating over light quarks" << endl;
  //cout << "evaluating tth hypo" << endl;
  ////integrate over the light quark angles
  //res = integrand->run( FinalState::LH, Hypothesis::TTH,  {PSVar::cos_q1, PSVar::phi_q1, PSVar::cos_qbar1, PSVar::phi_qbar1} );
  //cout << "p = " << res.p << " +- " << res.p_err << endl;
  //p0 = res.p;
  //
  //cout << "evaluating ttbb hypo" << endl;
  //res = integrand->run( FinalState::LH, Hypothesis::TTBB, {PSVar::cos_q1, PSVar::phi_q1, PSVar::cos_qbar1, PSVar::phi_qbar1} );
  //cout << "p = " << res.p << " +- " << res.p_err << endl;
  //p1 = res.p;
  //
  //mem_w = p0 / (p0 + 0.02*p1);
  //cout << "mem 022 discriminator " << mem_w << endl;

  integrand->next_event();

  delete integrand;
}
