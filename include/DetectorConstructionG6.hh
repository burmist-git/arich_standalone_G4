#ifndef DetectorConstructionG6_H
#define DetectorConstructionG6_H 1

//My
#include "DetectorConstruction.hh"

class DetectorConstructionG6 : public DetectorConstruction
{
public:
  
  DetectorConstructionG6() : DetectorConstruction()
  {
  }
     
  G4VPhysicalVolume* Construct();

};

#endif
