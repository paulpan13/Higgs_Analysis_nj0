////TMVA Classification macro for H->μμ analysis
//This macro implements the signal(Higgs) vs backgournd(Z boson)
//discrimination for the zero-jet channel (nj=0), following the ATLAS analysis methodology.
//MLP_vs_BDT_vs_DNN. The two methods will be tested and evaluated for their performance.
//Vol.2 Adding DNN.


#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

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

Int_t HiggsAnalysis0jet_Advanced(TString myMethodList = "") {}