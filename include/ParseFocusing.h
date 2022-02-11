
#ifndef PARSEFOCUSING_H
#define PARSEFOCUSING_H

#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/dk2nu.h"
#include "TrkChainVol.h"

#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>

namespace NeutrinoFocAna{
  
 /*! \class ParseFocusing
   *  \brief Classify in focusing components 
   *
   * This class classify the neutrino parent 
   * by the way that they has been focused:
   * underfocused, overfocused, H1 only, H2 only, 
   * unfocused and others.
   *
   */
  class ParseFocusing{
  public:
    
    //! default constructor
    ParseFocusing();
    //! Constructor using dk2nu
    void CalculateComponents(bsim::Dk2Nu* nu, bsim::DkMeta* meta, std::string beamline);
    //! Destructor
    ~ParseFocusing();
    std::string GetCategory();
    
    static ParseFocusing* getInstance();
    static void resetInstance();
    
  public:
    //NuMI
    //! Focusing name[6]: [Unfocused, H2only, H1only, Underfocused, Overfocused, Other]
    //LBNF
    //! Focusing name[6]: [??]
    std::string focname; 
  
  private:
    static ParseFocusing* instance;
  };
  
};

#endif 
