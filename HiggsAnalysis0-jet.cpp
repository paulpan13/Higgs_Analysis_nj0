//TMVA Classification macro for H->μμ analysis
//This macro implements the signal(Higgs) vs backgournd(Z boson)
//discrimination for the zero-jet channel (nj=0), following the ATLAS analysis methodology.
//MLP_vs_BDT. The two methods will be tested and evaluated for their performance.

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int HiggsAnalysis0jet(Tstring myMethodList = ""){
    // This initializes TMVA
   TMVA::Tools::Instance();
   //MVA methods
   std::map<std::string,int> Use;

   //Neural Networks
   Use["MLP"] = 1;

   //BDT
   Use["BDT"]=1;

   if (myMethodList!="") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return 1;
         }
         Use[regMethod] = 1;
      }
   }
   std::cout << "==> Start HiggsAnalysis0jet" << std::endl;

   TString outfileName("Higgs_0jet_Results.root");
   TFile* outputFile = TFile::Open(outfileName,"Recreate");

   //Factory & DataLoader
   TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G;D:AnalysisType=Classification");
   TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

   //pT
     dataloader->AddVariable("sqrt(lep_pt[0]*lep_pt[0] + lep_pt[1]*lep_pt[1] + 2*lep_pt[0]*lep_pt[1]*cos(lep_phi[0]-lep_phi[1]))", "pT_mumu", "MeV", 'F');

    //Rapidity (Y_mumu)
    dataloader->AddVariable("(lep_eta[0]+lep_eta[1])/2.0", "Y_mumu", "", 'F');

    // Γ. Γωνία Collins-Soper (cos_theta_star)
    //muon-pz (pz = pt * sinh(eta))
    TString pz0 = "(lep_pt[0]*sinh(lep_eta[0]))";
    TString pz1 = "(lep_pt[1]*sinh(lep_eta[1]))";
    TString pzLL = "(" + pz0 + "+" + pz1 + ")";
    
    // Light-cone components: P+ = E + pz, P- = E - pz
    TString Pp0 = "(lep_e[0]+" + pz0 + ")";
    TString Pm0 = "(lep_e[0]-" + pz0 + ")";
    TString Pp1 = "(lep_e[1]+" + pz1 + ")";
    TString Pm1 = "(lep_e[1]-" + pz1 + ")";
    
    //paper: 2*(P1+*P2- - P1-*P2+) / sqrt(m^2 * (m^2 + pT^2)) * sign(pz_ll)
    TString cosThetaFormula = "2*(" + Pp0 + "*" + Pm1 + " - " + Pm0 + "*" + Pp1 + ") / (sqrt(m_mumu*m_mumu * (m_mumu*m_mumu + pT_mumu*pT_mumu))) * (" + pzLL + "/abs(" + pzLL + "))";
    
    dataloader->AddVariable(cosThetaFormula, "cos_theta_star", "", 'F');

    //Spectators
    dataloader->AddSpectator("m_mumu", 'F'); 
    dataloader->AddSpectator("eventNumber", 'L');
    dataloader->AddSpectator("jet_n", 'I');

    //Signal(ggF,VBF) and Background(Z) files
    TFile *Signal_ggF1 = TFile::Open("ggF1.root");
    TFile *Signal_ggF2 = TFile::Open("ggF2.root");
    TFile *Signal_VBF = TFile::Open("VBF.root");
    TFile *Backg_Z = TFile::Open("Zskimmed_nj0_totalEntries.root");
    
    //Tree names = analysis for every File
    TTree *Tree_ggF1 = (TTree*) Signal_ggF1 -> Get("analysis");
    TTree *Tree_ggF2 = (TTree*) Signal_ggF2 -> Get("analysis");
    TTree *Tree_VBF = (TTree*) Signal_VBF -> Get("analysis");
    TTree *Tree_backgZ = (TTree*) Backg_Z -> Get("analysis");
    
    //We add them to dataloader and we set the weight for both Signal and Background to zero
    Double_t signal_weight = 1.0;
    Double_t backg_weight = 1.0;
   
    dataloader -> AddSignalTree(Tree_ggF1,signal_weight);
    dataloader -> AddSignalTree(Tree_ggF2,signal_weight);
    dataloader -> AddSignalTree(Tree_VBF,signal_weight);
    dataloader -> AddBackgroundTree(Tree_backgZ,backg_weight);
    
    //Cuts 
    TCut mycuts = "jet_n == 0 && m_mumu >= 120000 && m_mumu <=130000";
    TCut mycutb = mycuts;
    
   
    //Spliting data to Training and Testing
    dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

    //Booking Methods
    if(Use["MLP"])
    factory->BookMethod(dataloader,TMVA::Types::kMLP,"MLP","H:!V:NeuronType = tanh: VarTransform = N: Ncycles = 600: HiddenLayers = N+5: TestRate = 5 :!UseRegulator");

    //commment

}
