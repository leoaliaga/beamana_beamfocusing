
#include "TrkVol.h"

namespace NeutrinoFocAna{ 
  
  TrkVol::TrkVol(){
    
    for(int i=0; i<3; i++){
      P[i] = -999999;
      X[i] = -999999;
    }
    Pmag = -999999;
    Vol  = "NoDefinied";
    
  }
  
  TrkVol::TrkVol(double Mom[], double Pos[], std::string volname){
    
    for(int i=0; i<3; i++){
      P[i] = Mom[i]; 
      X[i] = Pos[i];
    }
    
    Pmag = std::sqrt(Mom[0]*Mom[0] + Mom[1]*Mom[1] +Mom[2]*Mom[2]);
    if( Pos[3] < -999998) Pmag = -999999.;
    
    //Volume:
    TrkVol::Vol = volname;
    
  }

  //----------------------------------------------------------------------
  TrkVol::~TrkVol(){
    
  }
  
  std::ostream& TrkVol::print(std::ostream& os) const {
    using namespace std;
    os<<"=> "<<setw(5)<<"p3 : ";
    for(int i=0; i<3; i++) {
      os<<setiosflags(ios::fixed) << setprecision(2)<<setw(6)<<P[i]<<" ";
    }
    os<<"=> "<<setw(5)<<"x3 :";
    for(int i=0; i<3; i++) {
      os<<setiosflags(ios::fixed) << setprecision(2)<<setw(6)<<X[i]<<" ";
    }
    os<<"=> Vol:"<<Vol<<endl;
    return os;
  }

}
