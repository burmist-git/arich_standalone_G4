#ifndef DetectorConstructionG5_H
#define DetectorConstructionG5_H 1

//My
#include "DetectorConstruction.hh"

class DetectorConstructionG5 : public DetectorConstruction
{
public:
  
  DetectorConstructionG5() : DetectorConstruction()
  {
  }
     
  G4VPhysicalVolume* Construct();

};

#endif
