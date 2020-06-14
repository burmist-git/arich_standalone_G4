//My
#include "DetectorConstructionG5.hh"
#include "SensitiveDetector.hh"
#include "MagneticField.hh"

//G4
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Paraboloid.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Color.hh"
#include "G4TwoVector.hh"
#include "G4SDManager.hh"
#include "globals.hh"
//magnetic field
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4UserLimits.hh"
//GDML
//#include <G4GDMLParser.hh>

//root 
#include "TMath.h"

#include "G4SystemOfUnits.hh"

G4VPhysicalVolume* DetectorConstructionG5::Construct()
{

  ConstructField();
  
  G4double world_sizeX = 300.0*cm;
  G4double world_sizeY = 300.0*cm;
  G4double world_sizeZ = 30.0*cm;

  // Box
  G4double volume1_sizeX  = 304.8*mm;
  G4double volume1_sizeY  =  95.3*mm;
  G4double volume1_sizeZ  =   5.0*mm;
  G4double volume1R       = 970.0*mm;
  G4double volume1Phi     = 11.25*deg;
  G4double volume1_X0     =  volume1R*cos(volume1Phi);
  G4double volume1_Y0     =  volume1R*sin(volume1Phi);
  G4double volume1_Z0     =  0.0*mm;

  // Tube to subtract
  G4double volume2_Rmin  = 0.0*mm;
  G4double volume2_Rmax  = 2.0*mm;


  // 
  // Define World Volume
  //
  G4VSolid *world_solid = new G4Box("World",world_sizeX/2.0,world_sizeY/2.0,world_sizeZ/2.0);
  G4LogicalVolume *world_logical = new G4LogicalVolume(world_solid,Air,"World");
  G4VPhysicalVolume *world_physical = new G4PVPlacement(0,G4ThreeVector(),world_logical,"World",0,false,0);

  //
  // volume1
  //
  G4VSolid *volume1_body_solid = new G4Box("volume1_body_solid", volume1_sizeX/2.0, volume1_sizeY/2.0, volume1_sizeZ/2.0);
  //G4LogicalVolume *volume1_body_logical = new G4LogicalVolume( volume1_body_solid, Air, "volume1_body_logical");

  //
  // volume2
  //
  G4VSolid *volume2_body_solid = new G4Tubs("volume2_body_solid", volume2_Rmin, volume2_Rmax, (volume1_sizeX+1)/2.0, 0, 360.0*deg);
  //G4LogicalVolume *volume2_body_logical = new G4LogicalVolume( volume2_body_solid, Air, "volume2_body_logical");


  G4RotationMatrix Ra_sub;
  G4ThreeVector Ta_sub;
  G4Transform3D Tr_sub;
  Ta_sub.setX(0.0);
  Ta_sub.setY(100.0/2.0-8);
  Ta_sub.setZ(3.0);
  Ra_sub.rotateY(90.0*deg);
  Tr_sub = G4Transform3D(Ra_sub, Ta_sub);
  G4SubtractionSolid* substraction_solid = new G4SubtractionSolid("substraction_solid", volume1_body_solid, volume2_body_solid, Tr_sub);
  for(int i = 1;i<6;i++){
    Ta_sub.setX(0.0);
    Ta_sub.setY(100.0/2.0 - 8 - 8*2*i);
    Ta_sub.setZ(3.0);
    Tr_sub = G4Transform3D(Ra_sub, Ta_sub);
    substraction_solid = new G4SubtractionSolid("substraction_solid", substraction_solid, volume2_body_solid, Tr_sub);
  }

  G4LogicalVolume *substraction_logical = new G4LogicalVolume( substraction_solid, Air, "substraction_logical");
  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;
  Ta.setX(volume1_X0);
  Ta.setY(volume1_Y0);
  Ta.setZ(volume1_Z0);
  Ra.rotateZ(volume1Phi);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *substraction_physical = new G4PVPlacement(Tr,                      //Transformation
							       substraction_logical,    //its logical volume
							       "substraction_physical", //its name
							       world_logical,           //its mother volume
							       false,                   //no boolean operation
							       0);                      //copy number

  /*
  for(int i = 0;i<2;i++){
    if (substraction){ 
      G4cout<<"i = "<<i<<G4endl;
      substraction 
    }
    else{
      G4cout<<"i = "<<i<<G4endl;
      substraction = new G4SubtractionSolid("Box+CylinderMoved", volume1_body_solid, volume2_body_solid, Tr);
    }
  }
  */

  return world_physical;
}
