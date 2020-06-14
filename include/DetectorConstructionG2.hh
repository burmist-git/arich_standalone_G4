#ifndef DetectorConstructionG2_H
#define DetectorConstructionG2_H 1

//My
#include "DetectorConstruction.hh"

class DetectorConstructionG2 : public DetectorConstruction
{
public:
  
  DetectorConstructionG2() : DetectorConstruction()
  {
  }
     
  G4VPhysicalVolume* Construct();

};

#endif
