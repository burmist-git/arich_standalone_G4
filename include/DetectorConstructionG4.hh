#ifndef DetectorConstructionG4_H
#define DetectorConstructionG4_H 1

//My
#include "DetectorConstruction.hh"

class DetectorConstructionG4 : public DetectorConstruction
{
public:
  
  DetectorConstructionG4() : DetectorConstruction()
  {
  }
     
  G4VPhysicalVolume* Construct();

};

#endif
