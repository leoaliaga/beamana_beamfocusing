#ifndef TRKCHAINVOL_H
#define TRKCHAINVOL_H

#include "TrkVol.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"

#include <vector>
#include <iostream>
#include <string>
#include <iomanip>

namespace NeutrinoFocAna{
  
  /*! \class TrkChainVol
   *  \brief Information about the chain of neutrino parent in volumes.
  */

  class TrkChainVol{
  public:

    //! boring old default constructor
    TrkChainVol();
    
    //! create an trajectory chain from the new dk2nu format
    TrkChainVol(bsim::Dk2Nu* nu, bsim::DkMeta* meta, std::string beamline);

    //! vector of neutrino parent trajectories
    std::vector<TrkVol> trkvol_chain;

    std::ostream& print(std::ostream& os=std::cout) const;
    
  private:
    
  };
  
  
}
#endif
