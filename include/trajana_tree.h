
#ifndef trajana_tree_h
#define trajana_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <string>

/*! Class to load the trajana_tree tree */
class trajana_tree{

 public: 

/** the constructor needs a TChain*/
  trajana_tree(TChain* chain);
  virtual ~trajana_tree();
  
  /** Get ntuple entry*/
  void GetEntry(Int_t ientry); 
  
  int    nupdg;
  int    parpdg;
  double  impwgt; 
  std::vector<double>*       enu     = 0;
  std::vector<double>*       geowgt  = 0;
  std::vector<std::string>*  prevol  = 0;
  std::vector<std::string>*  postvol = 0;
  std::vector<double>*       posx    = 0;
  std::vector<double>*       posy    = 0;
  std::vector<double>*       posz    = 0;
  std::vector<double>*       momx    = 0;
  std::vector<double>*       momy    = 0;
  std::vector<double>*       momz    = 0;

  int ntrees;
  int nentries;
  
 private:
  
  TChain* fChain;

};

#endif
