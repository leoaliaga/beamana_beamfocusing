
#include "trajana_tree.h"

trajana_tree::trajana_tree(TChain* chain){

  trajana_tree::fChain = chain;

  fChain->SetBranchAddress("nupdg"  , &nupdg);
  fChain->SetBranchAddress("parpdg" , &parpdg);
  fChain->SetBranchAddress("impwgt" , &impwgt);
  fChain->SetBranchAddress("enu"    , &enu);
  fChain->SetBranchAddress("geowgt" , &geowgt);
  fChain->SetBranchAddress("prevol" , &prevol );
  fChain->SetBranchAddress("postvol", &postvol);
  fChain->SetBranchAddress("posx"    ,&posx);
  fChain->SetBranchAddress("posy"    ,&posy);
  fChain->SetBranchAddress("posz"    ,&posz);
  fChain->SetBranchAddress("momx"    ,&momx);
  fChain->SetBranchAddress("momy"    ,&momy);
  fChain->SetBranchAddress("momz"    ,&momz);

  fChain->SetMakeClass(1);
  
  trajana_tree::ntrees   = fChain->GetNtrees();
  trajana_tree::nentries = fChain->GetEntries();
  
}

void trajana_tree::GetEntry(Int_t ientry){

  fChain->GetEntry(ientry);
  
}

trajana_tree::~trajana_tree(){
  
}
 
