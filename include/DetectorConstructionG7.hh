#ifndef DetectorConstructionG7_H
#define DetectorConstructionG7_H 1

//My
#include "DetectorConstruction.hh"

class DetectorConstructionG7 : public DetectorConstruction
{
public:
  
  DetectorConstructionG7() : DetectorConstruction()
  {
  }
     
  G4VPhysicalVolume* Construct();

};

#endif
