#include <algorithm>
#include <utility>
#include <cassert>
#include <cmath>

#include <Math/VectorUtil.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include "TTH/CommonClassifier/interface/DLBDTVars.h"




DLBDTMvaVariablesBase::DLBDTMvaVariablesBase():
eventWeight_(DLBDTMvaVariableFloat(name_eventWeight_))
{}



DLBDTMvaVariablesBase::DLBDTMvaVariablesBase(const double& eventWeight):
eventWeight_(DLBDTMvaVariableFloat(name_eventWeight_))
{
    // Event weight
    eventWeight_.setValue(eventWeight);
}



void DLBDTMvaVariablesBase::clearVariables(std::vector<DLBDTMvaVariablesBase*>& v_DLBDTMvaVariables){
    for(auto& DLBDTMvaVariables : v_DLBDTMvaVariables) delete DLBDTMvaVariables;
    v_DLBDTMvaVariables.clear();
}





//------------------------------------





DLBDTMvaVariablesEventClassification::DLBDTMvaVariablesEventClassification():
DLBDTMvaVariablesBase(),
multiplicity_jets_(DLBDTMvaVariableInt(name_multiplicity_jets_)),
btagDiscriminatorAverage_tagged_(DLBDTMvaVariableFloat(name_btagDiscriminatorAverage_tagged_)),
btagDiscriminatorAverage_untagged_(DLBDTMvaVariableFloat(name_btagDiscriminatorAverage_untagged_)),
minDeltaR_jet_jet_(DLBDTMvaVariableFloat(name_minDeltaR_jet_jet_)),
minDeltaR_tag_tag_(DLBDTMvaVariableFloat(name_minDeltaR_tag_tag_)),
avgDeltaR_jet_jet_(DLBDTMvaVariableFloat(name_avgDeltaR_jet_jet_)),
avgDeltaR_jet_tag_(DLBDTMvaVariableFloat(name_avgDeltaR_jet_tag_)),
avgDeltaR_tag_tag_(DLBDTMvaVariableFloat(name_avgDeltaR_tag_tag_)),
ptSum_jets_leptons_(DLBDTMvaVariableFloat(name_ptSum_jets_leptons_)),
multiplicity_higgsLikeDijet15_(DLBDTMvaVariableInt(name_multiplicity_higgsLikeDijet15_)),
mass_higgsLikeDijet_(DLBDTMvaVariableFloat(name_mass_higgsLikeDijet_)),
mass_higgsLikeDijet2_(DLBDTMvaVariableFloat(name_mass_higgsLikeDijet2_)),
mass_jet_jet_min_deltaR_(DLBDTMvaVariableFloat(name_mass_jet_jet_min_deltaR_)),
mass_tag_tag_min_deltaR_(DLBDTMvaVariableFloat(name_mass_tag_tag_min_deltaR_)),
mass_jet_tag_min_deltaR_(DLBDTMvaVariableFloat(name_mass_jet_tag_min_deltaR_)),
mass_tag_tag_max_mass_(DLBDTMvaVariableFloat(name_mass_tag_tag_max_mass_)),
median_mass_jet_jet_(DLBDTMvaVariableFloat(name_median_mass_jet_jet_)),
maxDeltaEta_jet_jet_(DLBDTMvaVariableFloat(name_maxDeltaEta_jet_jet_)),
maxDeltaEta_tag_tag_(DLBDTMvaVariableFloat(name_maxDeltaEta_tag_tag_)),
HT_jets_(DLBDTMvaVariableFloat(name_HT_jets_)),
HT_tags_(DLBDTMvaVariableFloat(name_HT_tags_)),
pT_jet_jet_min_deltaR_(DLBDTMvaVariableFloat(name_pT_jet_jet_min_deltaR_)),
pT_jet_tag_min_deltaR_(DLBDTMvaVariableFloat(name_pT_jet_tag_min_deltaR_)),
pT_tag_tag_min_deltaR_(DLBDTMvaVariableFloat(name_pT_tag_tag_min_deltaR_)),
mass_jet_jet_jet_max_pT_(DLBDTMvaVariableFloat(name_mass_jet_jet_jet_max_pT_)),
mass_jet_tag_tag_max_pT_(DLBDTMvaVariableFloat(name_mass_jet_tag_tag_max_pT_)),
centrality_jets_leps_(DLBDTMvaVariableFloat(name_centrality_jets_leps_)),
centrality_tags_(DLBDTMvaVariableFloat(name_centrality_tags_)),
twist_jet_jet_max_mass_(DLBDTMvaVariableFloat(name_twist_jet_jet_max_mass_)),
twist_jet_tag_max_mass_(DLBDTMvaVariableFloat(name_twist_jet_tag_max_mass_)),
twist_tag_tag_max_mass_(DLBDTMvaVariableFloat(name_twist_tag_tag_max_mass_)),
twist_tag_tag_min_deltaR_(DLBDTMvaVariableFloat(name_twist_tag_tag_min_deltaR_)),
sphericity_jet_(DLBDTMvaVariableFloat(name_sphericity_jet_)),
aplanarity_jet_(DLBDTMvaVariableFloat(name_aplanarity_jet_)),
circularity_jet_(DLBDTMvaVariableFloat(name_circularity_jet_)),
isotropy_jet_(DLBDTMvaVariableFloat(name_isotropy_jet_)),
C_jet_(DLBDTMvaVariableFloat(name_C_jet_)),
D_jet_(DLBDTMvaVariableFloat(name_D_jet_)),
transSphericity_jet_(DLBDTMvaVariableFloat(name_transSphericity_jet_)),
sphericity_tag_(DLBDTMvaVariableFloat(name_sphericity_tag_)),
aplanarity_tag_(DLBDTMvaVariableFloat(name_aplanarity_tag_)),
circularity_tag_(DLBDTMvaVariableFloat(name_circularity_tag_)),
isotropy_tag_(DLBDTMvaVariableFloat(name_isotropy_tag_)),
C_tag_(DLBDTMvaVariableFloat(name_C_tag_)),
D_tag_(DLBDTMvaVariableFloat(name_D_tag_)),
transSphericity_tag_(DLBDTMvaVariableFloat(name_transSphericity_tag_)),
H0_jet_(DLBDTMvaVariableFloat(name_H0_jet_)),
H1_jet_(DLBDTMvaVariableFloat(name_H1_jet_)),
H2_jet_(DLBDTMvaVariableFloat(name_H2_jet_)),
H3_jet_(DLBDTMvaVariableFloat(name_H3_jet_)),
H4_jet_(DLBDTMvaVariableFloat(name_H4_jet_)),
R1_jet_(DLBDTMvaVariableFloat(name_R1_jet_)),
R2_jet_(DLBDTMvaVariableFloat(name_R2_jet_)),
R3_jet_(DLBDTMvaVariableFloat(name_R3_jet_)),
R4_jet_(DLBDTMvaVariableFloat(name_R4_jet_)),
H0_tag_(DLBDTMvaVariableFloat(name_H0_tag_)),
H1_tag_(DLBDTMvaVariableFloat(name_H1_tag_)),
H2_tag_(DLBDTMvaVariableFloat(name_H2_tag_)),
H3_tag_(DLBDTMvaVariableFloat(name_H3_tag_)),
H4_tag_(DLBDTMvaVariableFloat(name_H4_tag_)),
R1_tag_(DLBDTMvaVariableFloat(name_R1_tag_)),
R2_tag_(DLBDTMvaVariableFloat(name_R2_tag_)),
R3_tag_(DLBDTMvaVariableFloat(name_R3_tag_)),
R4_tag_(DLBDTMvaVariableFloat(name_R4_tag_))
{}



DLBDTMvaVariablesEventClassification::DLBDTMvaVariablesEventClassification(
    const int multiplicity_jets,
    const double& btagDiscriminatorAverage_tagged,
    const double& btagDiscriminatorAverage_untagged,
    const double& minDeltaR_jet_jet,
    const double& minDeltaR_tag_tag,
    const double& avgDeltaR_jet_jet,
    const double& avgDeltaR_jet_tag,
    const double& avgDeltaR_tag_tag,
    const double& ptSum_jets_leptons,
    const int multiplicity_higgsLikeDijet15,
    const double& mass_higgsLikeDijet,
    const double& mass_higgsLikeDijet2,
    const double& mass_jet_jet_min_deltaR,
    const double& mass_tag_tag_min_deltaR,
    const double& mass_jet_tag_min_deltaR,
    const double& mass_tag_tag_max_mass,
    const double& median_mass_jet_jet,
    const double& maxDeltaEta_jet_jet,
    const double& maxDeltaEta_tag_tag,
    const double& HT_jets,
    const double& HT_tags,
    const double& pT_jet_jet_min_deltaR,
    const double& pT_jet_tag_min_deltaR,
    const double& pT_tag_tag_min_deltaR,
    const double& mass_jet_jet_jet_max_pT,
    const double& mass_jet_tag_tag_max_pT,
    const double& centrality_jets_leps,
    const double& centrality_tags,
    const double& twist_jet_jet_max_mass,
    const double& twist_jet_tag_max_mass,
    const double& twist_tag_tag_max_mass,
    const double& twist_tag_tag_min_deltaR,
    const double& sphericity_jet,
    const double& aplanarity_jet,
    const double& circularity_jet,
    const double& isotropy_jet,
    const double& C_jet,
    const double& D_jet,
    const double& transSphericity_jet,
    const double& sphericity_tag,
    const double& aplanarity_tag,
    const double& circularity_tag,
    const double& isotropy_tag,
    const double& C_tag,
    const double& D_tag,
    const double& transSphericity_tag,
    const double& H0_jet,
    const double& H1_jet,
    const double& H2_jet,
    const double& H3_jet,
    const double& H4_jet,
    const double& R1_jet,
    const double& R2_jet,
    const double& R3_jet,
    const double& R4_jet,
    const double& H0_tag,
    const double& H1_tag,
    const double& H2_tag,
    const double& H3_tag,
    const double& H4_tag,
    const double& R1_tag,
    const double& R2_tag,
    const double& R3_tag,
    const double& R4_tag):
DLBDTMvaVariablesBase(),
multiplicity_jets_(DLBDTMvaVariableInt(name_multiplicity_jets_)),
btagDiscriminatorAverage_tagged_(DLBDTMvaVariableFloat(name_btagDiscriminatorAverage_tagged_)),
btagDiscriminatorAverage_untagged_(DLBDTMvaVariableFloat(name_btagDiscriminatorAverage_untagged_)),
minDeltaR_jet_jet_(DLBDTMvaVariableFloat(name_minDeltaR_jet_jet_)),
minDeltaR_tag_tag_(DLBDTMvaVariableFloat(name_minDeltaR_tag_tag_)),
avgDeltaR_jet_jet_(DLBDTMvaVariableFloat(name_avgDeltaR_jet_jet_)),
avgDeltaR_jet_tag_(DLBDTMvaVariableFloat(name_avgDeltaR_jet_tag_)),
avgDeltaR_tag_tag_(DLBDTMvaVariableFloat(name_avgDeltaR_tag_tag_)),
ptSum_jets_leptons_(DLBDTMvaVariableFloat(name_ptSum_jets_leptons_)),
multiplicity_higgsLikeDijet15_(DLBDTMvaVariableInt(name_multiplicity_higgsLikeDijet15_)),
mass_higgsLikeDijet_(DLBDTMvaVariableFloat(name_mass_higgsLikeDijet_)),
mass_higgsLikeDijet2_(DLBDTMvaVariableFloat(name_mass_higgsLikeDijet2_)),
mass_jet_jet_min_deltaR_(DLBDTMvaVariableFloat(name_mass_jet_jet_min_deltaR_)),
mass_tag_tag_min_deltaR_(DLBDTMvaVariableFloat(name_mass_tag_tag_min_deltaR_)),
mass_jet_tag_min_deltaR_(DLBDTMvaVariableFloat(name_mass_jet_tag_min_deltaR_)),
mass_tag_tag_max_mass_(DLBDTMvaVariableFloat(name_mass_tag_tag_max_mass_)),
median_mass_jet_jet_(DLBDTMvaVariableFloat(name_median_mass_jet_jet_)),
maxDeltaEta_jet_jet_(DLBDTMvaVariableFloat(name_maxDeltaEta_jet_jet_)),
maxDeltaEta_tag_tag_(DLBDTMvaVariableFloat(name_maxDeltaEta_tag_tag_)),
HT_jets_(DLBDTMvaVariableFloat(name_HT_jets_)),
HT_tags_(DLBDTMvaVariableFloat(name_HT_tags_)),
pT_jet_jet_min_deltaR_(DLBDTMvaVariableFloat(name_pT_jet_jet_min_deltaR_)),
pT_jet_tag_min_deltaR_(DLBDTMvaVariableFloat(name_pT_jet_tag_min_deltaR_)),
pT_tag_tag_min_deltaR_(DLBDTMvaVariableFloat(name_pT_tag_tag_min_deltaR_)),
mass_jet_jet_jet_max_pT_(DLBDTMvaVariableFloat(name_mass_jet_jet_jet_max_pT_)),
mass_jet_tag_tag_max_pT_(DLBDTMvaVariableFloat(name_mass_jet_tag_tag_max_pT_)),
centrality_jets_leps_(DLBDTMvaVariableFloat(name_centrality_jets_leps_)),
centrality_tags_(DLBDTMvaVariableFloat(name_centrality_tags_)),
twist_jet_jet_max_mass_(DLBDTMvaVariableFloat(name_twist_jet_jet_max_mass_)),
twist_jet_tag_max_mass_(DLBDTMvaVariableFloat(name_twist_jet_tag_max_mass_)),
twist_tag_tag_max_mass_(DLBDTMvaVariableFloat(name_twist_tag_tag_max_mass_)),
twist_tag_tag_min_deltaR_(DLBDTMvaVariableFloat(name_twist_tag_tag_min_deltaR_)),
sphericity_jet_(DLBDTMvaVariableFloat(name_sphericity_jet_)),
aplanarity_jet_(DLBDTMvaVariableFloat(name_aplanarity_jet_)),
circularity_jet_(DLBDTMvaVariableFloat(name_circularity_jet_)),
isotropy_jet_(DLBDTMvaVariableFloat(name_isotropy_jet_)),
C_jet_(DLBDTMvaVariableFloat(name_C_jet_)),
D_jet_(DLBDTMvaVariableFloat(name_D_jet_)),
transSphericity_jet_(DLBDTMvaVariableFloat(name_transSphericity_jet_)),
sphericity_tag_(DLBDTMvaVariableFloat(name_sphericity_tag_)),
aplanarity_tag_(DLBDTMvaVariableFloat(name_aplanarity_tag_)),
circularity_tag_(DLBDTMvaVariableFloat(name_circularity_tag_)),
isotropy_tag_(DLBDTMvaVariableFloat(name_isotropy_tag_)),
C_tag_(DLBDTMvaVariableFloat(name_C_tag_)),
D_tag_(DLBDTMvaVariableFloat(name_D_tag_)),
transSphericity_tag_(DLBDTMvaVariableFloat(name_transSphericity_tag_)),
H0_jet_(DLBDTMvaVariableFloat(name_H0_jet_)),
H1_jet_(DLBDTMvaVariableFloat(name_H1_jet_)),
H2_jet_(DLBDTMvaVariableFloat(name_H2_jet_)),
H3_jet_(DLBDTMvaVariableFloat(name_H3_jet_)),
H4_jet_(DLBDTMvaVariableFloat(name_H4_jet_)),
R1_jet_(DLBDTMvaVariableFloat(name_R1_jet_)),
R2_jet_(DLBDTMvaVariableFloat(name_R2_jet_)),
R3_jet_(DLBDTMvaVariableFloat(name_R3_jet_)),
R4_jet_(DLBDTMvaVariableFloat(name_R4_jet_)),
H0_tag_(DLBDTMvaVariableFloat(name_H0_tag_)),
H1_tag_(DLBDTMvaVariableFloat(name_H1_tag_)),
H2_tag_(DLBDTMvaVariableFloat(name_H2_tag_)),
H3_tag_(DLBDTMvaVariableFloat(name_H3_tag_)),
H4_tag_(DLBDTMvaVariableFloat(name_H4_tag_)),
R1_tag_(DLBDTMvaVariableFloat(name_R1_tag_)),
R2_tag_(DLBDTMvaVariableFloat(name_R2_tag_)),
R3_tag_(DLBDTMvaVariableFloat(name_R3_tag_)),
R4_tag_(DLBDTMvaVariableFloat(name_R4_tag_))
{
    multiplicity_jets_.setValue(multiplicity_jets);
    btagDiscriminatorAverage_tagged_.setValue(btagDiscriminatorAverage_tagged);
    btagDiscriminatorAverage_untagged_.setValue(btagDiscriminatorAverage_untagged);
    minDeltaR_jet_jet_.setValue(minDeltaR_jet_jet);
    minDeltaR_tag_tag_.setValue(minDeltaR_tag_tag);
    avgDeltaR_jet_jet_.setValue(avgDeltaR_jet_jet);
    avgDeltaR_jet_tag_.setValue(avgDeltaR_jet_tag);
    avgDeltaR_tag_tag_.setValue(avgDeltaR_tag_tag);
    ptSum_jets_leptons_.setValue(ptSum_jets_leptons);
    multiplicity_higgsLikeDijet15_.setValue(multiplicity_higgsLikeDijet15);
    mass_higgsLikeDijet_.setValue(mass_higgsLikeDijet);
    mass_higgsLikeDijet2_.setValue(mass_higgsLikeDijet2);
    mass_jet_jet_min_deltaR_.setValue(mass_jet_jet_min_deltaR);
    mass_jet_tag_min_deltaR_.setValue(mass_jet_tag_min_deltaR);
    mass_tag_tag_min_deltaR_.setValue(mass_tag_tag_min_deltaR);
    mass_tag_tag_max_mass_.setValue(mass_tag_tag_max_mass);
    median_mass_jet_jet_.setValue(median_mass_jet_jet);
    maxDeltaEta_jet_jet_.setValue(maxDeltaEta_jet_jet);
    maxDeltaEta_tag_tag_.setValue(maxDeltaEta_tag_tag);
    HT_jets_.setValue(HT_jets);
    HT_tags_.setValue(HT_tags);
    pT_jet_jet_min_deltaR_.setValue(pT_jet_jet_min_deltaR);
    pT_jet_tag_min_deltaR_.setValue(pT_jet_tag_min_deltaR);
    pT_tag_tag_min_deltaR_.setValue(pT_tag_tag_min_deltaR);
    mass_jet_jet_jet_max_pT_.setValue(mass_jet_jet_jet_max_pT);
    mass_jet_tag_tag_max_pT_.setValue(mass_jet_tag_tag_max_pT);
    centrality_jets_leps_.setValue(centrality_jets_leps);
    centrality_tags_.setValue(centrality_tags);
    twist_jet_jet_max_mass_.setValue(twist_jet_jet_max_mass);
    twist_jet_tag_max_mass_.setValue(twist_jet_tag_max_mass);
    twist_tag_tag_max_mass_.setValue(twist_tag_tag_max_mass);
    twist_tag_tag_min_deltaR_.setValue(twist_tag_tag_min_deltaR);
    
    sphericity_jet_.setValue(sphericity_jet);
    aplanarity_jet_.setValue(aplanarity_jet);
    circularity_jet_.setValue(circularity_jet);
    isotropy_jet_.setValue(isotropy_jet);
    C_jet_.setValue(C_jet);
    D_jet_.setValue(D_jet);
    transSphericity_jet_.setValue(transSphericity_jet);
    
    sphericity_tag_.setValue(sphericity_tag);
    aplanarity_tag_.setValue(aplanarity_tag);
    circularity_tag_.setValue(circularity_tag);
    isotropy_tag_.setValue(isotropy_tag);
    C_tag_.setValue(C_tag);
    D_tag_.setValue(D_tag);
    transSphericity_tag_.setValue(transSphericity_tag);
    
    H0_jet_.setValue(H0_jet);
    H1_jet_.setValue(H1_jet);
    H2_jet_.setValue(H2_jet);
    H3_jet_.setValue(H3_jet);
    H4_jet_.setValue(H4_jet);
    R1_jet_.setValue(R1_jet);
    R2_jet_.setValue(R2_jet);
    R3_jet_.setValue(R3_jet);
    R4_jet_.setValue(R4_jet);
    
    H0_tag_.setValue(H0_tag);
    H1_tag_.setValue(H1_tag);
    H2_tag_.setValue(H2_tag);
    H3_tag_.setValue(H3_tag);
    H4_tag_.setValue(H4_tag);
    R1_tag_.setValue(R1_tag);
    R2_tag_.setValue(R2_tag);
    R3_tag_.setValue(R3_tag);
    R4_tag_.setValue(R4_tag);
}



DLBDTMvaVariablesEventClassification* DLBDTMvaVariablesEventClassification::fillVariables(
    const std::vector<TLorentzVector>& leptons, const std::vector<double>& selectedLeptonCharge,
    const std::vector<TLorentzVector>& jets, const std::vector<double>& jetBtags,
    const double btagWP)
{
    using ROOT::Math::VectorUtil::DeltaR;
    using ROOT::Math::VectorUtil::DeltaPhi;
    
    //build vectors of indices to be compatible with the rest if the code
    std::vector<int> jetIndices;
    std::vector<int> bjetIndices;
    for(unsigned int i = 0; i < jets.size(); ++i){
        jetIndices.push_back(i);
        if(jetBtags.at(i) >= btagWP) bjetIndices.push_back(i);
    }
    
    // Calculate lepton based quantities
    double sumPt_lepton(0.);
    double sumE_lepton(0.);
    for(const TLorentzVector& lepton : leptons){
        sumPt_lepton += lepton.Pt();
        sumE_lepton += lepton.E();
    }
    
    // Identify indices of jet pairs fulfilling given criterion
    std::pair<int, int> p_minDeltaR_jet_jet;
    std::pair<int, int> p_minDeltaR_jet_tag;
    std::pair<int, int> p_minDeltaR_tag_tag;
    std::pair<int, int> p_maxDeltaEta_jet_jet;
    std::pair<int, int> p_maxDeltaEta_jet_tag;
    std::pair<int, int> p_maxDeltaEta_tag_tag;
    std::pair<int, int> p_maxMass_jet_jet;
    std::pair<int, int> p_maxMass_jet_tag;
    std::pair<int, int> p_maxMass_tag_tag;
    std::pair<int, int> p_closestHiggsMass_jet_jet;
    std::pair<int, int> p_closestHiggsMass_jet_tag;
    std::pair<int, int> p_closestHiggsMass_tag_tag;
    
    // Values of jet pairs fulfilling given criterion
    double minDeltaR_jet_jet(999.);
    double minDeltaR_jet_tag(999.);
    double minDeltaR_tag_tag(999.);
    double maxDeltaEta_jet_jet(-999.);
    double maxDeltaEta_jet_tag(-999.);
    double maxDeltaEta_tag_tag(-999.);
    double maxMass_jet_jet(-999.);
    double maxMass_jet_tag(-999.);
    double maxMass_tag_tag(-999.);
    double closestHiggsMass_jet_jet(-999.);
    double closestHiggsMass_jet_tag(-999.);
    double closestHiggsMass_tag_tag(-999.);
    
    int multiplicity_jet_jet(0);
    int multiplicity_jet_tag(0);
    int multiplicity_tag_tag(0);
    int multiplicity_higgsLike15_jet_tag(0);
    double sumDeltaR_jet_jet(0.);
    double sumDeltaR_jet_tag(0.);
    double sumDeltaR_tag_tag(0.);
    
    // Sums of per-jet quantities
    double sumBtagDiscriminator_tagged(0.);
    double sumBtagDiscriminator_untagged(0.);
    double sumPt_jet(0.);
    double sumE_jet(0.);
    double sumPt_tag(0.);
    double sumE_tag(0.);
    
    // Loop over all jets and jet pairs
    for(auto i_index = jetIndices.begin(); i_index != jetIndices.end(); ++i_index){
        const bool btagged1 = std::find(bjetIndices.begin(), bjetIndices.end(), *i_index) != bjetIndices.end();
        const TLorentzVector& jet1 = jets.at(*i_index);
        
        sumPt_jet += jet1.Pt();
        sumE_jet += jet1.E();
        
        const double& btagDiscriminator(jetBtags.at(*i_index));
        // Set values<0. to 0.
        const double btagDiscriminatorPositive(btagDiscriminator>=0. ? btagDiscriminator : 0.);
        // Avoid b-tag values where the algorithm did not work, giving values>1., setting them to 1.
        const double btagDiscriminatorInRange(btagDiscriminatorPositive<=1. ? btagDiscriminatorPositive : 1.);
        if(btagged1){
            sumBtagDiscriminator_tagged += btagDiscriminatorInRange;
            sumPt_tag += jet1.Pt();
            sumE_tag += jet1.E();
        }
        else{
            sumBtagDiscriminator_untagged += btagDiscriminatorInRange;
        }
        
        for(auto j_index = i_index + 1; j_index != jetIndices.end(); ++j_index){
            const bool btagged2 = std::find(bjetIndices.begin(),bjetIndices.end(), *j_index) != bjetIndices.end();
            const TLorentzVector& jet2 = jets.at(*j_index);
            
            constexpr double higgsMass(125.);
            const double deltaR = DeltaR(jet1, jet2);
            const double deltaEta = std::abs(jet1.Eta() - jet2.Eta());
            const double mass = (jet1 + jet2).M();
            
            
            // Jet-jet combinations
            if(deltaR < minDeltaR_jet_jet){
                minDeltaR_jet_jet = deltaR;
                p_minDeltaR_jet_jet = std::make_pair(*i_index, *j_index);
            }
            if(deltaEta > maxDeltaEta_jet_jet){
                maxDeltaEta_jet_jet = deltaEta;
                p_maxDeltaEta_jet_jet = std::make_pair(*i_index, *j_index);
            }
            if(mass > maxMass_jet_jet){
                maxMass_jet_jet = mass;
                p_maxMass_jet_jet = std::make_pair(*i_index, *j_index);
            }
            if(std::abs(higgsMass - mass) < std::abs(higgsMass - closestHiggsMass_jet_jet)){
                closestHiggsMass_jet_jet = mass;
                p_closestHiggsMass_jet_jet = std::make_pair(*i_index, *j_index);
            }
            
            ++multiplicity_jet_jet;
            sumDeltaR_jet_jet += deltaR;
            
            
            // Jet-tag combinations
            if(!btagged1 && !btagged2) continue;
            if(deltaR < minDeltaR_jet_tag){
                minDeltaR_jet_tag = deltaR;
                p_minDeltaR_jet_tag = std::make_pair(*i_index, *j_index);
            }
            if(deltaEta > maxDeltaEta_jet_tag){
                maxDeltaEta_jet_tag = deltaEta;
                p_maxDeltaEta_jet_tag = std::make_pair(*i_index, *j_index);
            }
            if(mass > maxMass_jet_tag){
                maxMass_jet_tag = mass;
                p_maxMass_jet_tag = std::make_pair(*i_index, *j_index);
            }
            if(std::abs(higgsMass - mass) < std::abs(higgsMass - closestHiggsMass_jet_tag)){
                closestHiggsMass_jet_tag = mass;
                p_closestHiggsMass_jet_tag = std::make_pair(*i_index, *j_index);
            }
            
            ++multiplicity_jet_tag;
            if(std::abs(higgsMass - mass) < 15.) ++multiplicity_higgsLike15_jet_tag;
            sumDeltaR_jet_tag += deltaR;
            
            
            // Tag-tag combinations
            if(!btagged1 || !btagged2) continue;
            if(deltaR < minDeltaR_tag_tag){
                minDeltaR_tag_tag = deltaR;
                p_minDeltaR_tag_tag = std::make_pair(*i_index, *j_index);
            }
            if(deltaEta > maxDeltaEta_tag_tag){
                maxDeltaEta_tag_tag = deltaEta;
                p_maxDeltaEta_tag_tag = std::make_pair(*i_index, *j_index);
            }
            if(mass > maxMass_tag_tag){
                maxMass_tag_tag = mass;
                p_maxMass_tag_tag = std::make_pair(*i_index, *j_index);
            }
            if(std::abs(higgsMass - mass) < std::abs(higgsMass - closestHiggsMass_tag_tag)){
                closestHiggsMass_tag_tag = mass;
                p_closestHiggsMass_tag_tag = std::make_pair(*i_index, *j_index);
            }
            
            ++multiplicity_tag_tag;
            sumDeltaR_tag_tag += deltaR;
        }
    }
    
    
    // Calculate jet quantities
    const int multiplicity_jet(jetIndices.size());
    const int multiplicity_tag(bjetIndices.size());
    const int multiplicity_untagged(multiplicity_jet - multiplicity_tag);
    const double averageBtagDiscriminator_tagged = multiplicity_tag>0 ? sumBtagDiscriminator_tagged/static_cast<double>(multiplicity_tag) : 0.;
    const double averageBtagDiscriminator_untagged = multiplicity_untagged>0 ? sumBtagDiscriminator_untagged/static_cast<double>(multiplicity_untagged) : 0.;
    const double sumPt_jet_lepton = sumPt_lepton + sumPt_jet;
    const double centrality_jet_lepton = (sumPt_lepton + sumPt_jet)/(sumE_lepton + sumE_jet);
    const double centrality_btag = sumPt_tag/sumE_tag;
    
    
    // Calculate jet pair quantities
    const double averageDeltaR_jet_jet = multiplicity_jet_jet>0 ? sumDeltaR_jet_jet/multiplicity_jet_jet : 0.;
    const double averageDeltaR_jet_tag = multiplicity_jet_tag>0 ? sumDeltaR_jet_tag/multiplicity_jet_tag : 0.;
    const double averageDeltaR_tag_tag = multiplicity_tag_tag>0 ? sumDeltaR_tag_tag/multiplicity_tag_tag : 0.;
    const double pt_minDeltaR_jet_jet = (jets.at(p_minDeltaR_jet_jet.first) + jets.at(p_minDeltaR_jet_jet.second)).Pt();
    const double pt_minDeltaR_jet_tag = (jets.at(p_minDeltaR_jet_tag.first) + jets.at(p_minDeltaR_jet_tag.second)).Pt();
    const double pt_minDeltaR_tag_tag = (jets.at(p_minDeltaR_tag_tag.first) + jets.at(p_minDeltaR_tag_tag.second)).Pt();
    const double mass_minDeltaR_jet_jet = (jets.at(p_minDeltaR_jet_jet.first) + jets.at(p_minDeltaR_jet_jet.second)).M();
    const double mass_minDeltaR_jet_tag = (jets.at(p_minDeltaR_jet_tag.first) + jets.at(p_minDeltaR_jet_tag.second)).M();
    const double mass_minDeltaR_tag_tag = (jets.at(p_minDeltaR_tag_tag.first) + jets.at(p_minDeltaR_tag_tag.second)).M();
    const double mass_maxMass_tag_tag = (jets.at(p_maxMass_tag_tag.first) + jets.at(p_maxMass_tag_tag.second)).M();
    const double twist_maxMass_jet_jet = TMath::ATan(DeltaPhi(jets.at(p_maxMass_jet_jet.first), jets.at(p_maxMass_jet_jet.second))/(jets.at(p_maxMass_jet_jet.first).Eta() - jets.at(p_maxMass_jet_jet.second).Eta()));
    const double twist_maxMass_jet_tag = TMath::ATan(DeltaPhi(jets.at(p_maxMass_jet_tag.first), jets.at(p_maxMass_jet_tag.second))/(jets.at(p_maxMass_jet_tag.first).Eta() - jets.at(p_maxMass_jet_tag.second).Eta()));
    const double twist_maxMass_tag_tag = TMath::ATan(DeltaPhi(jets.at(p_maxMass_tag_tag.first), jets.at(p_maxMass_tag_tag.second))/(jets.at(p_maxMass_tag_tag.first).Eta() - jets.at(p_maxMass_tag_tag.second).Eta()));
    const double twist_minDeltaR_tag_tag = TMath::ATan(DeltaPhi(jets.at(p_minDeltaR_tag_tag.first), jets.at(p_minDeltaR_tag_tag.second))/(jets.at(p_minDeltaR_tag_tag.first).Eta() - jets.at(p_minDeltaR_tag_tag.second).Eta()));
    
    
    // Calculate median mass
    std::vector<double> v_mass_jet_jet;
    for(auto i_index = jetIndices.begin(); i_index != jetIndices.end(); ++i_index){
        const TLorentzVector& jet1 = jets.at(*i_index);
        for(auto j_index = i_index + 1; j_index != jetIndices.end(); ++j_index){
            const TLorentzVector& jet2 = jets.at(*j_index);
            const double mass = (jet1 + jet2).M();
            v_mass_jet_jet.push_back(mass);
        }
    }
    std::sort(v_mass_jet_jet.begin(), v_mass_jet_jet.end());
    
    double medianMass_jet_jet(-999.);
    if(v_mass_jet_jet.size()){
        const size_t nJetJet(v_mass_jet_jet.size());
        if(nJetJet%2 == 0) medianMass_jet_jet = (v_mass_jet_jet.at(nJetJet/2) + v_mass_jet_jet.at(nJetJet/2-1))/2;
        else medianMass_jet_jet = v_mass_jet_jet.at((nJetJet-1)/2);
    }
    
    
    // Calculate tri-jet masses
    double mass_maxPt_jet_jet_jet(-999.);
    double mass_maxPt_jet_tag_tag(-999.);
    double maxPt_jet_jet_jet(-999.);
    double maxPt_jet_tag_tag(-999.);
    for(auto i_index = jetIndices.begin(); i_index != jetIndices.end(); ++i_index){
        const TLorentzVector& jet1  = jets.at(*i_index);
        const bool tagged1 = std::find(bjetIndices.begin(),bjetIndices.end(), *i_index) != bjetIndices.end();
        for(auto j_index = i_index + 1; j_index != jetIndices.end(); ++j_index){
            const TLorentzVector& jet2  = jets.at(*j_index);
            const bool tagged2 = std::find(bjetIndices.begin(),bjetIndices.end(), *j_index) != bjetIndices.end();
            for(auto k_index = j_index + 1; k_index != jetIndices.end(); ++k_index){
                const TLorentzVector& jet3  = jets.at(*k_index);
                const bool tagged3 = std::find(bjetIndices.begin(),bjetIndices.end(), *k_index) != bjetIndices.end();
                const double pt_jet_jet_jet = (jet1+jet2+jet3).Pt();
                if(pt_jet_jet_jet > maxPt_jet_jet_jet){
                    maxPt_jet_jet_jet = pt_jet_jet_jet;
                    mass_maxPt_jet_jet_jet = (jet1+jet2+jet3).M();
                }
                if((!tagged1 && !tagged2) || (!tagged1 && !tagged3) || (!tagged2 && !tagged3)) continue;
                if(pt_jet_jet_jet > maxPt_jet_tag_tag){
                    maxPt_jet_tag_tag = pt_jet_jet_jet;
                    mass_maxPt_jet_tag_tag = (jet1+jet2+jet3).M();
                }
            }
        }
    }
    
    
    // Event shape variable for jets
    const EventShapeVariables eventshape_jets(jets);
    
    const double sphericity_jet = eventshape_jets.sphericity();
    const double aplanarity_jet = eventshape_jets.aplanarity();
    const double circularity_jet = eventshape_jets.circularity();
    const double isotropy_jet = eventshape_jets.isotropy();
    const double C_jet = eventshape_jets.C();
    const double D_jet = eventshape_jets.D();
    const double transSphericity_jet = eventshape_jets.transSphericity();
    
    const double H0_jet = eventshape_jets.H(0);
    const double H1_jet = eventshape_jets.H(1);
    const double H2_jet = eventshape_jets.H(2);
    const double H3_jet = eventshape_jets.H(3);
    const double H4_jet = eventshape_jets.H(4);
    
    const double R1_jet = eventshape_jets.R(1);
    const double R2_jet = eventshape_jets.R(2);
    const double R3_jet = eventshape_jets.R(3);
    const double R4_jet = eventshape_jets.R(4);
    
    
    // Event shape variables for b-tagged jets
    std::vector<TLorentzVector> btags;
    for(const int index : bjetIndices) btags.push_back(jets.at(index));
    
    const EventShapeVariables eventshape_tags(btags);        
    
    const double sphericity_tag = eventshape_tags.sphericity();
    const double aplanarity_tag = eventshape_tags.aplanarity();
    const double circularity_tag = eventshape_tags.circularity();
    const double isotropy_tag = eventshape_tags.isotropy();
    const double C_tag = eventshape_tags.C();
    const double D_tag = eventshape_tags.D();
    const double transSphericity_tag = eventshape_tags.transSphericity();
    
    const double H0_tag = eventshape_tags.H(0);
    const double H1_tag = eventshape_tags.H(1);
    const double H2_tag = eventshape_tags.H(2);
    const double H3_tag = eventshape_tags.H(3);
    const double H4_tag = eventshape_tags.H(4);
    
    const double R1_tag = eventshape_tags.R(1);
    const double R2_tag = eventshape_tags.R(2);
    const double R3_tag = eventshape_tags.R(3);
    const double R4_tag = eventshape_tags.R(4);
    
    
    // Fill variables
    return new DLBDTMvaVariablesEventClassification(multiplicity_jet,
                                                    averageBtagDiscriminator_tagged,
                                                    averageBtagDiscriminator_untagged,
                                                    minDeltaR_jet_jet,
                                                    minDeltaR_tag_tag,
                                                    averageDeltaR_jet_jet,
                                                    averageDeltaR_jet_tag,
                                                    averageDeltaR_tag_tag,
                                                    sumPt_jet_lepton,
                                                    multiplicity_higgsLike15_jet_tag,
                                                    closestHiggsMass_jet_jet,
                                                    closestHiggsMass_jet_tag,
                                                    mass_minDeltaR_jet_jet,
                                                    mass_minDeltaR_tag_tag,
                                                    mass_minDeltaR_jet_tag,
                                                    mass_maxMass_tag_tag,
                                                    medianMass_jet_jet,
                                                    maxDeltaEta_jet_jet,
                                                    maxDeltaEta_tag_tag,
                                                    sumPt_jet,
                                                    sumPt_tag,
                                                    pt_minDeltaR_jet_jet,
                                                    pt_minDeltaR_jet_tag,
                                                    pt_minDeltaR_tag_tag,
                                                    mass_maxPt_jet_jet_jet,
                                                    mass_maxPt_jet_tag_tag,
                                                    centrality_jet_lepton,
                                                    centrality_btag,
                                                    twist_maxMass_jet_jet,
                                                    twist_maxMass_jet_tag,
                                                    twist_maxMass_tag_tag,
                                                    twist_minDeltaR_tag_tag,
                                                    sphericity_jet,
                                                    aplanarity_jet,
                                                    circularity_jet,
                                                    isotropy_jet,
                                                    C_jet,
                                                    D_jet,
                                                    transSphericity_jet,
                                                    sphericity_tag,
                                                    aplanarity_tag,
                                                    circularity_tag,
                                                    isotropy_tag,
                                                    C_tag,
                                                    D_tag,
                                                    transSphericity_tag,
                                                    H0_jet,
                                                    H1_jet,
                                                    H2_jet,
                                                    H3_jet,
                                                    H4_jet,
                                                    R1_jet,
                                                    R2_jet,
                                                    R3_jet,
                                                    R4_jet,
                                                    H0_tag,
                                                    H1_tag,
                                                    H2_tag,
                                                    H3_tag,
                                                    H4_tag,
                                                    R1_tag,
                                                    R2_tag,
                                                    R3_tag,
                                                    R4_tag);
}







// ---------------------------------- Class DLBDTMvaVariablesEventClassification::EventShapeVariables -------------------------------------------



DLBDTMvaVariablesEventClassification::EventShapeVariables::EventShapeVariables(const std::vector<TLorentzVector>& inputVectors):
inputVectors_(this->transformInputVectors(inputVectors))
{}



std::vector<ROOT::Math::XYZVector> DLBDTMvaVariablesEventClassification::EventShapeVariables::transformInputVectors(const std::vector<TLorentzVector>& jets)const
{
    std::vector<ROOT::Math::XYZVector> result;
    for(size_t i = 0; i < jets.size(); ++i){
        const TLorentzVector& jetP4(jets.at(i));
        const ROOT::Math::XYZVector jetP3 = ROOT::Math::XYZVector(jetP4.Px(), jetP4.Py(), jetP4.Pz());
        result.push_back(jetP3);
    }
    return result;
}



double DLBDTMvaVariablesEventClassification::EventShapeVariables::isotropy(unsigned int numberOfSteps)const
{
    double eIn(-1.);
    double eOut(-1.);
    
    const double deltaPhi = 2*TMath::Pi()/numberOfSteps;
    double phi(0.);
    for(unsigned int i = 0; i < numberOfSteps; ++i){
        phi += deltaPhi;
        
        // Sum over inner product of unit vectors and momenta
        double sum(0.);
        for(unsigned int j = 0; j < inputVectors_.size(); ++j)
            sum += TMath::Abs(TMath::Cos(phi)*inputVectors_.at(j).x()+TMath::Sin(phi)*inputVectors_.at(j).y());
        
        if(eOut<0. || sum<eOut) eOut = sum;
        if(eIn<0. || sum>eIn) eIn = sum;
    }
    
    return (eIn-eOut)/eIn;
}



double DLBDTMvaVariablesEventClassification::EventShapeVariables::circularity(unsigned int numberOfSteps)const
{
    double circularity(-1.);
    
    double area(0.);
    for(unsigned int i = 0; i < inputVectors_.size(); ++i)
        area += TMath::Sqrt(inputVectors_.at(i).x()*inputVectors_.at(i).x()+inputVectors_.at(i).y()*inputVectors_.at(i).y());
    
    const double deltaPhi = 2*TMath::Pi()/numberOfSteps;
    double phi(0.);
    for(unsigned int i = 0; i < numberOfSteps; ++i){
        phi += deltaPhi;
        
        double sum(0.);
        for(unsigned int j = 0; j < inputVectors_.size(); ++j)
            sum += TMath::Abs(TMath::Cos(phi)*inputVectors_.at(j).x() + TMath::Sin(phi)*inputVectors_.at(j).y());
        
        const double tmpCircularity = TMath::Pi()/2.*sum/area;
        if(circularity<0. || tmpCircularity<circularity) circularity = tmpCircularity;
    }
    
    return circularity;
}



double DLBDTMvaVariablesEventClassification::EventShapeVariables::sphericity(const double& r)const
{
    const TVectorD eigenValues = this->compEigenValues(r);
    return 1.5*(eigenValues(1) + eigenValues(2));
}



double DLBDTMvaVariablesEventClassification::EventShapeVariables::aplanarity(const double& r)const
{
    const TVectorD eigenValues = this->compEigenValues(r);
    return 1.5*eigenValues(2);
}



double DLBDTMvaVariablesEventClassification::EventShapeVariables::C(const double& r)const
{
    const TVectorD eigenValues = this->compEigenValues(r);
    return 3.*(eigenValues(0)*eigenValues(1) + eigenValues(0)*eigenValues(2) + eigenValues(1)*eigenValues(2));
}



double DLBDTMvaVariablesEventClassification::EventShapeVariables::D(const double& r)const
{
    const TVectorD eigenValues = this->compEigenValues(r);
    return 27.*eigenValues(0)*eigenValues(1)*eigenValues(2);
}



double DLBDTMvaVariablesEventClassification::EventShapeVariables::transSphericity(const double& r)const
{
    const TVectorD eigenValues = this->compEigenValues(r);
    return 2.*eigenValues(1)/(eigenValues(0) + eigenValues(1));
}



double DLBDTMvaVariablesEventClassification::EventShapeVariables::H(const int order)const
{
    const FoxWolframMoments fwm(inputVectors_, order);
    return fwm.H(order);
}



double DLBDTMvaVariablesEventClassification::EventShapeVariables::R(const int order)const
{
    const FoxWolframMoments fwm(inputVectors_, order);
    return fwm.R(order);
}



TVectorD DLBDTMvaVariablesEventClassification::EventShapeVariables::compEigenValues(const double& r)const
{
    TVectorD eigenValues(3);
    const TMatrixDSym tensor = this->compMomentumTensor(r);
    if(tensor.IsSymmetric() && tensor.NonZeros()!=0) tensor.EigenVectors(eigenValues);
    return eigenValues;
}



TMatrixDSym DLBDTMvaVariablesEventClassification::EventShapeVariables::compMomentumTensor(const double& r)const
{
    TMatrixDSym momentumTensor(3);
    momentumTensor.Zero();
    if(inputVectors_.size() < 2) return momentumTensor;
    
    // Fill momentumTensor from inputVectors
    double norm(1.);
    for(int i = 0; i < (int)inputVectors_.size(); ++i){
        const ROOT::Math::XYZVector& inputVector(inputVectors_.at(i));
        const double p2 = inputVector.Dot(inputVector);
        const double pR = r==2. ? p2 : TMath::Power(p2, 0.5*r);
        norm += pR;
        const double pRminus2 = r==2. ? 1. : TMath::Power(p2, 0.5*r - 1.);
        momentumTensor(0,0) += pRminus2*inputVector.x()*inputVector.x();
        momentumTensor(0,1) += pRminus2*inputVector.x()*inputVector.y();
        momentumTensor(0,2) += pRminus2*inputVector.x()*inputVector.z();
        momentumTensor(1,0) += pRminus2*inputVector.y()*inputVector.x();
        momentumTensor(1,1) += pRminus2*inputVector.y()*inputVector.y();
        momentumTensor(1,2) += pRminus2*inputVector.y()*inputVector.z();
        momentumTensor(2,0) += pRminus2*inputVector.z()*inputVector.x();
        momentumTensor(2,1) += pRminus2*inputVector.z()*inputVector.y();
        momentumTensor(2,2) += pRminus2*inputVector.z()*inputVector.z();
    }
    
    return (1./norm)*momentumTensor;
}







// ---------------------------------- Class DLBDTMvaVariablesEventClassification::FoxWolframMoments -------------------------------------------

DLBDTMvaVariablesEventClassification::FoxWolframMoments::FoxWolframMoments(const std::vector<ROOT::Math::XYZVector>& inputVectors,
                                                                           const int maxorder):
nMoment_(maxorder+1),
fwArray_(nMoment_),
sumArray_(nMoment_)
{
    for(int i = 0; i < nMoment_; ++i){ 
        fwArray_(i) = 0.;
        sumArray_(i) = 0.;
    }
    
    this->compute(inputVectors);
}



void DLBDTMvaVariablesEventClassification::FoxWolframMoments::compute(const std::vector<ROOT::Math::XYZVector>& inputVectors)
{
    if(inputVectors.size() == 0) return;
    
    // Loop over the all particle candidates
    double s(0.);
    for(size_t i = 0; i<inputVectors.size(); ++i){
        // Candidate particle's 3-momentum
        const TVector3 p1(inputVectors.at(i).X(), inputVectors.at(i).Y(), inputVectors.at(i).Z());
        const double pmag1 = p1.Mag();

        // Loop over other particle's candidates, starting at the next one in the list
        for(size_t j = i; j < inputVectors.size(); ++j){
            // Candidate particle's 3-momentum
            const TVector3 p2(inputVectors.at(j).X(),inputVectors.at(j).Y(),inputVectors.at(j).Z());
            
            // Cosine of the angle between the two candidate particles
            const double cosPhi = TMath::Cos(p1.Angle(p2));
            
            // Contribution of this pair of track
            // (note the factor 2 : the pair enters the sum twice)
            const double pmag2 = p2.Mag();
            for(int l = 0; l < nMoment_; ++l) sumArray_(l) += 2 * pmag1 * pmag2 * this->legendre(l, 0, cosPhi);
        }
        
        // Contribution for this moment
        for(int l = 0; l < nMoment_; ++l) sumArray_(l) += pmag1*pmag1*this->legendre(l, 0, 1.);
        
        // Total energy
        s += p1.Mag();
    }

    if(s <= 0.) return;

    // Normalize Fox Wolfram Moments
    for(int i = 0; i < nMoment_; ++i) fwArray_(i) = sumArray_(i)/TMath::Power(s, 2);
}



double DLBDTMvaVariablesEventClassification::FoxWolframMoments::R(const int order)const
{
    if(this->H(0)>0. && order<nMoment_) return (this->H(order)/this->H(0));
    return 0.;
}



double DLBDTMvaVariablesEventClassification::FoxWolframMoments::legendre(const int l, const int m, const double& x)const
{
    assert(m >= 0.);
    assert(m <= l);
    assert(std::abs(x) <= 1.);
    
    double pmm(1.);
    if(m > 0){
        const double somx2 = TMath::Sqrt((1. - x) * (1. + x));
        double fact = 1.;
        for(int i = 0; i < m; ++i){
            pmm *= -fact * somx2;
            fact += 2.;
        }
    }
    if(l == m) return pmm;

    double pmmp1 = x * (2 * m + 1) * pmm;
    if(l == m + 1) return pmmp1;
    
    for(int ll = m+2; ll <= l; ++ll){
        const double pll = (x*(2 * ll - 1)*pmmp1 - (ll + m - 1)*pmm)/(ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pmmp1;
}




