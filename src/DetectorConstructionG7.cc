//My
#include "DetectorConstructionG7.hh"
#include "SensitiveDetector.hh"
#include "MagneticField.hh"
#include "libxmlarichdata.h"
#include "tessellatedSolidStr.h"

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
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "globals.hh"
//magnetic field
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4UserLimits.hh"
//GDML
//#include <G4GDMLParser.hh>
//VGM
#include "Geant4GM/volumes/Factory.h"
#include "RootGM/volumes/Factory.h"
#include "TGeoManager.h"

//root 
#include "TMath.h"

#include "G4SystemOfUnits.hh"

G4VPhysicalVolume* DetectorConstructionG7::Construct()
{
  
  tessellatedSolidStr tessellatedSolid = xmlarichdata::readTessellatedSolidVerticesToDATfile("tessellatedSolidVertices.dat");
  //tessellatedSolid.printInfo();
 
  //xmlarichdata::quartz_Babar_DIRC_polishing_quality_bounds();

  ConstructField();
  
  G4double world_sizeX = 50.0*cm;
  G4double world_sizeY = 50.0*cm;
  G4double world_sizeZ = 50.0*cm;

  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;

  // 
  // Define World Volume
  //
  G4VSolid *world_solid = new G4Box("World",world_sizeX/2.0,world_sizeY/2.0,world_sizeZ/2.0);
  G4LogicalVolume *world_logical = new G4LogicalVolume(world_solid,Air,"World");
  G4VPhysicalVolume *world_physical = new G4PVPlacement(0,G4ThreeVector(),world_logical,"World",0,false,0);


  
  //
  //Tessellated volume
  //
  G4TessellatedSolid *volume_solid = new G4TessellatedSolid("TessellatedSolid");

  G4ThreeVector point_1;
  G4ThreeVector point_2;
  G4ThreeVector point_3;

  for(unsigned int i=0; i < tessellatedSolid.nCells; i++){
    //
    point_1.setX(tessellatedSolid.posV1[0][i]);
    point_1.setY(tessellatedSolid.posV1[1][i]);
    point_1.setZ(tessellatedSolid.posV1[2][i]);
    //
    point_2.setX(tessellatedSolid.posV2[0][i]);
    point_2.setY(tessellatedSolid.posV2[1][i]);
    point_2.setZ(tessellatedSolid.posV2[2][i]);
    //
    point_3.setX(tessellatedSolid.posV3[0][i]);
    point_3.setY(tessellatedSolid.posV3[1][i]);
    point_3.setZ(tessellatedSolid.posV3[2][i]);    
    //
    G4TriangularFacet * facet = new G4TriangularFacet(point_1, point_2, point_3, ABSOLUTE);
    volume_solid->AddFacet((G4VFacet*) facet);
    
  }

  volume_solid->SetSolidClosed(true);

  G4LogicalVolume *volume_logical = new G4LogicalVolume( volume_solid, Aluminum, "volume_logical");
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *volume_physical = new G4PVPlacement(Tr,                //Transformation
							 volume_logical,    //its logical volume                            
							 "volume_physical", //its name
							 world_logical,     //its mother  volume
							 false,             //no boolean operation
							 0);                //copy number
  
  

  //
  //Simple Box
  //
  G4double box_size = 5.0 * cm;
  G4VSolid *box_solid = new G4Box("box_solid", box_size/2.0, box_size/2.0, box_size/2.0);
  G4LogicalVolume *box_logical = new G4LogicalVolume( box_solid, Aluminum, "box_logical");
  Ta.setX(10 * cm);
  Ta.setY(10 * cm);
  Ta.setZ(0);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *box_physical = new G4PVPlacement(Tr,             //Transformation
						      box_logical,    //its logical volume                            
						      "box_physical", //its name
						      world_logical,  //its mother  volume
						      false,          //no boolean operation
						      0);             //copy number

  
  G4double maxStep   = 0.1*mm;
  G4double maxLength = 2.0*m;
  G4double maxTime   = 20.0*ns; 
  G4double minEkin   = 1.0/100*MeV;
  G4double mionRang  = 0.01*mm;
  stepLimit = new G4UserLimits(maxStep,maxLength,maxTime,minEkin,mionRang);
  
  // ---------------------------------------------------------------------------
  // VGM demo 
  // Export geometry in Root and save it in a file
  // 
  // 
  // Import Geant4 geometry to VGM
  Geant4GM::Factory g4Factory;
  g4Factory.SetIgnore(1);
  g4Factory.SetDebug(1);
  g4Factory.Import(world_physical);
  // 
  // Export VGM geometry to Root
  RootGM::Factory rtFactory;
  rtFactory.SetDebug(1);
  rtFactory.SetIgnore(1);
  g4Factory.Export(&rtFactory);
  gGeoManager->CloseGeometry();
  gGeoManager->Export("geometry.root");  
  //
  // end VGM demo
  //---------------------------------------------------------------------------

  return world_physical;
}
