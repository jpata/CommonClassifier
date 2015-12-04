#ifndef MEMCLASSIFIER_H
#define MEMCLASSIFIER_H

class MEMResult {
public:
  //likelihood ratio
  double p;

  //individual signal and background probabilities
  double p_sig;
  double p_bkg;

  //Integration uncertainties of the probabilities
  double p_err_sig;
  double p_err_bkg;
  
  //Integration uncertainties of the probabilities
  double n_perm_sig;
  double n_perm_bkg;
};

class MEMClassifier{
  
public:

  //The constructor loads the transfer functions and b-tag PDF-s
  MEMClassifier();
  ~MEMClassifier();

  // Call this method to return the BDT output, provide all necessary inputs. Jet CSV should be sorted the same way as jet p4. 
  // We could also write a class to contain the jet CSV and p4 information
  MEMResult GetOutput(
    const std::vector<TLorentzVector>& selectedLeptonP4,
    const std::vector<TLorentzVector>& selectedJetP4,
    const std::vector<double>& selectedJetCSV,
    const std::vector<TLorentzVector>& looseSelectedJetP4,
    const std::vector<double>& looseSelectedJetCSV,
    TLorentzVector& metP4
  );

  // returns the category of the last evaluated Event
  std::string GetCategoryOfLastEvaluation() const;

private:

  //Holds the transfer functions
  TFile* transfers;

  //This is the MEM workhorse from the ETH side
  Integrand* integrand; 
};

#endif
