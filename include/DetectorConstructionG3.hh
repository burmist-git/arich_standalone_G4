#ifndef DetectorConstructionG3_H
#define DetectorConstructionG3_H 1

//My
#include "DetectorConstruction.hh"

class DetectorConstructionG3 : public DetectorConstruction
{
public:
  
  DetectorConstructionG3() : DetectorConstruction()
  {
  }
     
  G4VPhysicalVolume* Construct();

};

#endif
