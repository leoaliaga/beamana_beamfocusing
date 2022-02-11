#ifndef TRKVOL_H
#define TRKVOL_H

#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace NeutrinoFocAna{
  
    /*! \class TrkVol
     *  \brief The information about a neutrino parent in a volume.   
 */
  class TrkVol{ 
 
  public:

    //!Default Constructor
    TrkVol();

    //! Constructor given kinematic of the particle.
    TrkVol(double Mom[], double Pos[], std::string volname);

    virtual ~TrkVol();
    
    //! Momentum vector of the particle
    double P[3]; 
    double Pmag;
    
    //! position vector of the particle
    double X[3]; 

    //! trajectory volume 
    std::string Vol;

    std::ostream& print(std::ostream& os) const;

  private:
    
  };

  
}
#endif
