#ifndef BDTCLASSIFIER_H
#define BDTCLASSIFIER_H

class CommonClassifier{
  
public:
  BDTClassifier();
  ~BDTClassifier();


  // Call this method to return the BDT output, provide all necessary inputs. Jet CSV should be sorted the same way as jet p4. 
  // We could also write a class to contain the jet CSV and p4 information
  double GetBDTOutput(const std::vector<TLorentzVector>& selectedLeptonP4, const std::vector<TLorentzVector>& selectedJetP4, const std::vector<double>& selectedJetCSV, const std::vector<TLorentzVector>& looseSelectedJetP4, const std::vector<double>& looseSelectedJetCSV, TLorentzVector& metP4);
  
  // returns the category of the last evaluated Event
  std::string GetCategoryOfLastEvaluation() const;
  
  // return the variable names and their values for the last evaluated event
  std::map<std::string,double> GetVariablesOfLastEvaluation() const;
  
  
private:
   //map holding the various BDT input variables
   std::map<std::string,double> variableMap;
  
   std::string categoryOfLastEvaluation;
   
};

#endif
