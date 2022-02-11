
#include "TrkChainVol.h"

namespace NeutrinoFocAna{ 
  TrkChainVol::TrkChainVol(){}
  
  TrkChainVol::TrkChainVol(bsim::Dk2Nu* nu, bsim::DkMeta* meta, std::string beamline){
    double defval = -999998;
    
    if(beamline == "NuMI"){
      for(int ii=0;ii<=7;ii++){
	bool not_wanted = (nu->traj[ii].trkz < defval);
	if(not_wanted)continue;
	double this_pos[3] = {nu->traj[ii].trkx ,nu->traj[ii].trky ,nu->traj[ii].trkz};
	double this_mom[3] = {nu->traj[ii].trkpx,nu->traj[ii].trkpy,nu->traj[ii].trkpz};
	std::string volname = "noname";
	////
	if(ii==0) volname = "BornTgt";
	if(ii==1) volname = "ExitTgt";
	if(ii==2){
	  if(nu->traj[ii].trkz > abs(defval))volname = "Neck1";
	  else                           volname = "EnterHorn1";
	}
	if(ii==3) volname = "ExitHorn1";
	if(ii==4){
	  if(nu->traj[ii].trkz > abs(defval)) volname = "Neck2";
	  else                           volname = "EnterHorn2";
	}
	if(ii==5) volname = "ExitHorn2";
	if(ii==6) volname = "EnterDVOL";
	if(ii==7) volname = "Decay";
	TrkVol trkvol(this_mom,this_pos,volname);
	trkvol_chain.push_back(trkvol);	
      }
      
    } //end of NuMI
    
    else if(beamline == "LBNF"){
      std::cout<<"=> Not implemented yet"<<std::endl;
    }
    else{
      std::cout<<"=> Only working with (NuMI/LBNF). You entered: "<<beamline<<std::endl;
      exit (1);
    }
    
  }

  std::ostream& TrkChainVol::print(std::ostream& os) const{
    using namespace std;
    os<<"==== TrkChainVol ===="<<std::endl;
    os<<"\n *trk vols*\n";
    for(size_t i=0; i<trkvol_chain.size(); i++){
      os<<"   ";
      trkvol_chain[i].print(os);
    }
    os<<endl;
    return os;
  }
  
}
