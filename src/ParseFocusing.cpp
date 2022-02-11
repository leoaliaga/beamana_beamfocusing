#include "ParseFocusing.h"

namespace NeutrinoFocAna{

  ParseFocusing* ParseFocusing::instance = 0;
  
  ParseFocusing::ParseFocusing(){
    ParseFocusing::focname = "none";
  }
  ////////////////////
  void ParseFocusing::CalculateComponents(bsim::Dk2Nu* nu, bsim::DkMeta* meta, 
					  std::string beamline){
    
    TrkChainVol tcv(nu,meta, beamline);
    
    ParseFocusing::focname = "Other";
    
    double mval    = -999999;
    double auxH1   = -999999;
    double auxH2   = -999999;
    bool bornInTGT = false;
    bool neck1     = false;
    bool neck2     = false;
    bool horn1     = false;
    bool horn2     = false;

    for(unsigned int itrk = 0;itrk<tcv.trkvol_chain.size();itrk++){
      if(tcv.trkvol_chain[itrk].Vol == "BornTgt")   bornInTGT= true;
      if(tcv.trkvol_chain[itrk].Vol == "Neck1")     neck1    = true;
      if(tcv.trkvol_chain[itrk].Vol == "Neck2")     neck2    = true;
      if(tcv.trkvol_chain[itrk].Vol == "EnterHorn1"){
	horn1    = true;
	mval  = -1.0*( tcv.trkvol_chain[itrk].X[0] )/( tcv.trkvol_chain[itrk].X[1] );
	auxH1 = ( tcv.trkvol_chain[itrk].X[1] ) - ( tcv.trkvol_chain[itrk].X[0] )*mval;
      }
      if(tcv.trkvol_chain[itrk].Vol == "EnterHorn2"){
	horn2    = true;
     	auxH2 = (tcv.trkvol_chain[itrk].X[1]) - ( tcv.trkvol_chain[itrk].X[0] )*mval;
      }
    }

    if(bornInTGT){
      if(neck1 && neck2) ParseFocusing::focname = "Unfocused";
      if(neck1 && horn2) ParseFocusing::focname = "H2only";    
      if(horn1 && neck2) ParseFocusing::focname = "H1only";
      if(horn1 && horn2){
	bool no_cross = (auxH1>0. && auxH2>0.) || (auxH1<0. && auxH2<0.);
	if(no_cross)ParseFocusing::focname = "Underfocused";
	else	    ParseFocusing::focname = "Overfocused";
	}
      }
  } 

  std::string ParseFocusing::GetCategory(){
    return focname;
  }
  //////////////////
  ParseFocusing::~ParseFocusing(){ 
  }
  
  ParseFocusing* ParseFocusing::getInstance(){
    if (instance == 0) instance = new ParseFocusing;
    return instance;
  }
  ////
  void ParseFocusing::resetInstance()
  {
    delete instance; 
    instance = 0; 
  }
  
};
