//My
#include "DetectorConstruction.hh"
#include "SensitiveDetector.hh"
#include "MagneticField.hh"
#include "libxmlarichdata.h"

//G4
#include "G4SystemOfUnits.hh"
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

using namespace CLHEP;

DetectorConstruction::DetectorConstruction()
{
  magField = new MagneticField();
  worldVisAtt = new G4VisAttributes();
  quartzVisAtt = new G4VisAttributes();
  aerogelVisAtt = new G4VisAttributes();
  sensitiveVisAtt = new G4VisAttributes();
  pmtboxVisAtt = new G4VisAttributes();
  absVisAtt = new G4VisAttributes();
  // Define Materials to be used
  DefineMaterials();
}

DetectorConstruction::~DetectorConstruction()
{
  delete magField;
  delete worldVisAtt;
  delete quartzVisAtt;
  delete sensitiveVisAtt;
  delete pmtboxVisAtt;
  delete absVisAtt;
  delete stepLimit;
}

void DetectorConstruction::xmlarichdataUtil(){
  ///////////////////////////////
  //    ARICH AEROGEL MODULE   //
  const G4int numAerogel = 36;
  G4double wl_min = 240; //nm
  G4double wl_max = 700; //nm
  G4double pe_min = xmlarichdata::lambda_nm_2_ev(wl_max); //e.p.
  G4double pe_max = xmlarichdata::lambda_nm_2_ev(wl_min); //e.p.
  G4double PhotonEnergyAerogel[numAerogel];
  G4double PhotonLambda[numAerogel];
  G4double AerogelRefractiveIndexA_lam[numAerogel];
  G4double AerogelRefractiveIndexB_lam[numAerogel];
  G4double AerogelRefractiveIndexA[numAerogel];
  G4double AerogelRefractiveIndexB[numAerogel];
  G4double AerogelRayleighScattLA[numAerogel];
  G4double AerogelRayleighScattLB[numAerogel];
  //xmlarichdata::hapdQE();
  //xmlarichdata::refractiveIndexSilicaAerogelPrintTestData();
  for(int i = 0; i < numAerogel; i++){
    PhotonLambda[i] = wl_min + (wl_max - wl_min)/(numAerogel-1)*i;
    AerogelRefractiveIndexA_lam[i] = xmlarichdata::refractiveIndexSilicaAerogel("A",PhotonLambda[i]);
    AerogelRefractiveIndexB_lam[i] = xmlarichdata::refractiveIndexSilicaAerogel("B",PhotonLambda[i]);
  }
  TGraph *gr_AerogelRefractiveIndex_vs_lam_A = new TGraph(numAerogel,PhotonLambda,AerogelRefractiveIndexA_lam);
  TGraph *gr_AerogelRefractiveIndex_vs_lam_B = new TGraph(numAerogel,PhotonLambda,AerogelRefractiveIndexB_lam);
  TGraph *gr_AerogelRefractiveIndex_vs_pe_A = xmlarichdata::convert_TGraph_Lambda_vs_pe(gr_AerogelRefractiveIndex_vs_lam_A);
  TGraph *gr_AerogelRefractiveIndex_vs_pe_B = xmlarichdata::convert_TGraph_Lambda_vs_pe(gr_AerogelRefractiveIndex_vs_lam_B);
  //gr_AerogelRefractiveIndex_vs_lam_A->SaveAs("gr_AerogelRefractiveIndex_vs_lam_A.root");
  //gr_AerogelRefractiveIndex_vs_lam_B->SaveAs("gr_AerogelRefractiveIndex_vs_lam_B.root");
  for(int i = 0; i < numAerogel; i++){
    PhotonEnergyAerogel[i] = pe_min + (pe_max - pe_min)/(numAerogel-1)*i;
    AerogelRefractiveIndexA[i] = gr_AerogelRefractiveIndex_vs_pe_A->Eval(PhotonEnergyAerogel[i]);
    AerogelRefractiveIndexB[i] = gr_AerogelRefractiveIndex_vs_pe_B->Eval(PhotonEnergyAerogel[i]);
  }
  TGraph *gr_AerogelRefractiveIndex_vs_peNormalStep_A = new TGraph( numAerogel, PhotonEnergyAerogel, AerogelRefractiveIndexA);
  TGraph *gr_AerogelRefractiveIndex_vs_peNormalStep_B = new TGraph( numAerogel, PhotonEnergyAerogel, AerogelRefractiveIndexB);
  gr_AerogelRefractiveIndex_vs_pe_A->SaveAs("gr_AerogelRefractiveIndex_vs_pe_A.root");
  gr_AerogelRefractiveIndex_vs_pe_B->SaveAs("gr_AerogelRefractiveIndex_vs_pe_B.root");
  gr_AerogelRefractiveIndex_vs_peNormalStep_A->SaveAs("gr_AerogelRefractiveIndex_vs_peNormalStep_A.root");
  gr_AerogelRefractiveIndex_vs_peNormalStep_B->SaveAs("gr_AerogelRefractiveIndex_vs_peNormalStep_B.root");
  //xmlarichdata::print_TGraph_info(gr_AerogelRefractiveIndex_vs_pe_A, "g4");
  //xmlarichdata::print_TGraph_info(gr_AerogelRefractiveIndex_vs_pe_B, "g4");
  xmlarichdata::print_TGraph_info(gr_AerogelRefractiveIndex_vs_peNormalStep_A, "g4");
  xmlarichdata::print_TGraph_info(gr_AerogelRefractiveIndex_vs_peNormalStep_B, "g4");
  ///////////////////////////////

  TGraph *gr_Aerogel_Rayleigh_scat_L_vs_lam_A = xmlarichdata::plots("aerogel_transmittance_vs_photon_wavelength_aerogel_xml_ver3.root", "gr_A_L_norm" , "" , "!-v");
  TGraph *gr_Aerogel_Rayleigh_scat_L_vs_lam_B = xmlarichdata::plots("aerogel_transmittance_vs_photon_wavelength_aerogel_xml_ver3.root", "gr_B_L_norm" , "" , "!-v");
  gr_Aerogel_Rayleigh_scat_L_vs_lam_A->SaveAs("gr_Aerogel_Rayleigh_scat_L_vs_lam_A.root");
  gr_Aerogel_Rayleigh_scat_L_vs_lam_B->SaveAs("gr_Aerogel_Rayleigh_scat_L_vs_lam_B.root");
  TGraph *gr_Aerogel_Rayleigh_scat_L_vs_pe_A = xmlarichdata::convert_TGraph_Lambda_vs_pe(gr_Aerogel_Rayleigh_scat_L_vs_lam_A);
  TGraph *gr_Aerogel_Rayleigh_scat_L_vs_pe_B = xmlarichdata::convert_TGraph_Lambda_vs_pe(gr_Aerogel_Rayleigh_scat_L_vs_lam_B);
  gr_Aerogel_Rayleigh_scat_L_vs_pe_A->SaveAs("gr_Aerogel_Rayleigh_scat_L_vs_pe_A.root");
  gr_Aerogel_Rayleigh_scat_L_vs_pe_B->SaveAs("gr_Aerogel_Rayleigh_scat_L_vs_pe_B.root");
  for(int i = 0; i < numAerogel; i++){
    AerogelRayleighScattLA[i] = gr_Aerogel_Rayleigh_scat_L_vs_pe_A->Eval(PhotonEnergyAerogel[i]);
    AerogelRayleighScattLB[i] = gr_Aerogel_Rayleigh_scat_L_vs_pe_B->Eval(PhotonEnergyAerogel[i]);
  }
  TGraph *gr_AerogelRayleighScattL_vs_peNormalStep_A = new TGraph( numAerogel, PhotonEnergyAerogel, AerogelRayleighScattLA);
  TGraph *gr_AerogelRayleighScattL_vs_peNormalStep_B = new TGraph( numAerogel, PhotonEnergyAerogel, AerogelRayleighScattLB);
  gr_AerogelRayleighScattL_vs_peNormalStep_A->SaveAs("gr_AerogelRayleighScattL_vs_peNormalStep_A.root");
  gr_AerogelRayleighScattL_vs_peNormalStep_B->SaveAs("gr_AerogelRayleighScattL_vs_peNormalStep_B.root");
  G4cout<<"gr_AerogelRayleighScattL_vs_peNormalStep_A"<<G4endl;
  xmlarichdata::print_TGraph_info(gr_AerogelRayleighScattL_vs_peNormalStep_A, "g4");
  G4cout<<"gr_AerogelRayleighScattL_vs_peNormalStep_B"<<G4endl;
  xmlarichdata::print_TGraph_info(gr_AerogelRayleighScattL_vs_peNormalStep_B, "g4");

}

void DetectorConstruction::DefineMaterials()
{

  G4String symbol;
  G4double a, z, density;
  G4int ncomponents, natoms;
  G4double fractionmass;

  // Define elements
  G4Element* H = 
    new G4Element("Hydrogen", symbol = "H", z = 1., a = 1.01*g/mole);
  G4Element* C = 
    new G4Element("Carbon",   symbol = "C", z = 6., a = 12.01*g/mole);
  G4Element* N = 
    new G4Element("Nitrogen", symbol = "N", z = 7., a = 14.01*g/mole);
  G4Element* O =
    new G4Element("Oxygen",   symbol = "O", z = 8., a = 16.00*g/mole);
  G4Element* Si = 
    new G4Element("Silicon",  symbol = "Si", z = 14., a = 28.09*g/mole);
  G4Element* Al = 
    new G4Element("Aluminum", symbol = "Al", z = 13., a = 26.98*g/mole);
  G4Element* F = 
    new G4Element("Fluorine", symbol = "F", z = 9., a = 19.0*g/mole);

  // Quartz Material (SiO2_cladd)
  SiO2_cladd = new G4Material("quartzCladd", density = 2.200*g/cm3, ncomponents = 2);
  SiO2_cladd->AddElement(Si, natoms = 1);
  SiO2_cladd->AddElement(O , natoms = 2);

  // Quartz Material (SiO2_coat)
  SiO2_coat = new G4Material("quartzCoat", density = 2.200*g/cm3, ncomponents = 2);
  SiO2_coat->AddElement(Si, natoms = 1);
  SiO2_coat->AddElement(O , natoms = 2);

  // Quartz Material (SiO2)
  SiO2 = new G4Material("quartz", density = 2.200*g/cm3, ncomponents = 2);
  SiO2->AddElement(Si, natoms = 1);
  SiO2->AddElement(O , natoms = 2);

  // Aerogel Material (SiO2)
  Aerogel = new G4Material("Aerogel", density = 2*g/cm3, ncomponents = 2);
  Aerogel->AddElement(Si, natoms = 1);
  Aerogel->AddElement(O , natoms = 2);

  // AerogelA Material (SiO2)
  //<Element fraction="0.665">O</Element> 
  //<Element fraction="0.042">H</Element>
  //<Element fraction="0.292">Si</Element> 
  //<Element fraction="0.001">C</Element> 
  //</Components> 
  AerogelA = new G4Material("AerogelA", density = 0.163478*g/cm3, ncomponents = 4);
  AerogelA->AddElement(O , fractionmass = 0.665);
  AerogelA->AddElement(H , fractionmass = 0.042);
  AerogelA->AddElement(Si, fractionmass = 0.292);
  AerogelA->AddElement(C , fractionmass = 0.001);

  // AerogelB Material (SiO2)
  AerogelB = new G4Material("AerogelB", density = 0.195491*g/cm3, ncomponents = 4);
  AerogelB->AddElement(O , fractionmass = 0.665);
  AerogelB->AddElement(H , fractionmass = 0.042);
  AerogelB->AddElement(Si, fractionmass = 0.292);
  AerogelB->AddElement(C , fractionmass = 0.001);

  // C4F10
  C4F10 = new G4Material("fluorocarbon", density = 1.8/1000*g/cm3, ncomponents = 2);
  C4F10->AddElement(C, natoms = 4);
  C4F10->AddElement(F, natoms = 10);

  // Air
  Air = new G4Material("Air", density = 1.290*mg/cm3, ncomponents = 2);
  Air->AddElement(N, fractionmass = 0.7);
  Air->AddElement(O, fractionmass = 0.3);

  // Aluminum
  Aluminum = new G4Material("Aluminum", density = 2.7*g/cm3, ncomponents = 1);
  Aluminum->AddElement(Al, fractionmass = 1.0);

  // Aluminum of the mirror
  AluminumMirr = new G4Material("mirrAluminum", density = 2.7*g/cm3, ncomponents = 1);
  AluminumMirr->AddElement(Al, fractionmass = 1.0);

  /*
  // Assign Materials
  world.material = Air;
  secA.material = SiO2;
  secB.material = SiO2;
  secC.material = SiO2;
  secWin.material = SiO2;
  fiberCorr.material = SiO2;
  fiberClad.material = SiO2_cladd;
  //fiberCoat.material = SiO2_coat;
  fiberCoat.material = Aluminum;
  sensitive.material = Aluminum;
  abs1.material = Aluminum;
  */

  //
  // Generate and Add Material Properties Table
  //						
  const G4int num = 36;
  G4double WaveLength[num];
  G4double Absorption[num];      // Default value for absorption
  G4double AirAbsorption[num];
  G4double AirRefractiveIndex[num];
  G4double PhotonEnergy[num];

  // Absorption of quartz per 1m
  G4double QuartzAbsorption[num] =
    {0.999572036,0.999544661,0.999515062,0.999483019,0.999448285,
     0.999410586,0.999369611,0.999325013,0.999276402,0.999223336,
     0.999165317,0.999101778,0.999032079,0.998955488,0.998871172,
     0.998778177,0.99867541 ,0.998561611,0.998435332,0.998294892,
     0.998138345,0.997963425,0.997767484,0.997547418,0.99729958 ,
     0.99701966 ,0.99670255 ,0.996342167,0.995931242,0.995461041,
     0.994921022,0.994298396,0.993577567,0.992739402,0.991760297,
     0.990610945};

  G4double C4F10RefractiveIndex[num];
  G4double AerogelRefractiveIndex[num];

  for (int i=0; i<num; i++) {
    WaveLength[i] = (300 + i*10)*nanometer;
    // Aerogel and C4F10
    Absorption[i] = 100*m;
    // Air
    AirAbsorption[i] = 4000.*cm;
    AirRefractiveIndex[i] = 1.;
    // C4F10
    C4F10RefractiveIndex[i] = 1.0014;
    // Aerogel
    //AerogelRefractiveIndex[i] = 1.03;
    AerogelRefractiveIndex[i] = 1.055;
    PhotonEnergy[num - (i+1)] = twopi*hbarc/WaveLength[i];
    // Absorption is given per length and G4 needs mean free path
    // length, calculate it here
    // mean free path length - taken as probablility equal 1/e
    // that the photon will be absorbed
    QuartzAbsorption[i] = (-1)/log(QuartzAbsorption[i])*100*cm;
    //EpotekAbsorption[i] = (-1)/log(EpotekAbsorption[i])*
    //epotekBarJoint.thickness;
  }

  G4double QuartzRefractiveIndex[num] =
    {1.456535,1.456812,1.4571  ,1.457399,1.457712,1.458038,
     1.458378,1.458735,1.459108,1.4595  ,1.459911,1.460344,
     1.460799,1.46128 ,1.461789,1.462326,1.462897,1.463502,
     1.464146,1.464833,1.465566,1.46635 ,1.46719 ,1.468094,
     1.469066,1.470116,1.471252,1.472485,1.473826,1.475289,
     1.476891,1.478651,1.480592,1.482739,1.485127,1.487793};

  G4double CladdingRefractiveIndex[num];

  for(int i=0; i<num; i++){
    CladdingRefractiveIndex[i] = TMath::Sqrt(QuartzRefractiveIndex[i]*QuartzRefractiveIndex[i]-0.22*0.22); 
  }

  ///////////////////////////////
  //    ARICH AEROGEL MODULE   //
  const G4int numAerogel = 36;
  G4double PhotonEnergyAerogel[numAerogel] = 
    {1.77120*eV, 1.86820*eV, 1.96519*eV, 2.06219*eV, 2.15918*eV, 2.25618*eV,
     2.35317*eV, 2.45016*eV, 2.54716*eV, 2.64415*eV, 2.74115*eV, 2.83814*eV,
     2.93514*eV, 3.03213*eV, 3.12913*eV, 3.22612*eV, 3.32312*eV, 3.42011*eV,
     3.51710*eV, 3.61410*eV, 3.71109*eV, 3.80809*eV, 3.90508*eV, 4.00208*eV,
     4.09907*eV, 4.19607*eV, 4.29306*eV, 4.39005*eV, 4.48705*eV, 4.58404*eV,
     4.68104*eV, 4.77803*eV, 4.87503*eV, 4.97202*eV, 5.06902*eV, 5.16601*eV};
  G4double AerogelRefractiveIndexA[numAerogel] = 
    {1.04489, 1.04496, 1.04504, 1.04513, 1.04522, 1.04531,
     1.04541, 1.04551, 1.04562, 1.04573, 1.04584, 1.04597,
     1.04609, 1.04622, 1.04636, 1.04650, 1.04665, 1.04680,
     1.04696, 1.04712, 1.04729, 1.04747, 1.04764, 1.04783,
     1.04802, 1.04822, 1.04843, 1.04864, 1.04886, 1.04908,
     1.04931, 1.04956, 1.04980, 1.05006, 1.05032, 1.05058};
  G4double AerogelRefractiveIndexB[numAerogel] = 
    {1.05372, 1.05381, 1.05390, 1.05400, 1.05410, 1.05421,
     1.05433, 1.05445, 1.05457, 1.05470, 1.05484, 1.05498,
     1.05512, 1.05528, 1.05544, 1.05560, 1.05577, 1.05595,
     1.05613, 1.05632, 1.05652, 1.05672, 1.05693, 1.05714,
     1.05737, 1.05760, 1.05784, 1.05808, 1.05833, 1.05859,
     1.05886, 1.05914, 1.05942, 1.05972, 1.06002, 1.06033};
  G4double AerogelRayleighScattLA[numAerogel] =
    { 312.3800, 268.67500, 229.16200, 197.05200, 169.61700, 146.19500,
      126.4130, 109.41200,  94.91120,  82.47880,  71.98050,  63.22030,
       55.6823,  49.14980,  43.43180,  38.44740,  34.10230,  30.31140,
       26.9924,  24.03170,  21.38230,  19.05140,  16.99240,  15.12560,
       13.4597,  12.03510,  10.76940,   9.63790,   8.64043,   7.72909,
       7.01191,   6.32322,   5.73803,   5.16746,   4.69441,   4.22136};
  G4double AerogelRayleighScattLB[numAerogel] = 
    { 232.24900, 200.51600, 171.77900, 148.07500, 127.78400, 110.30400,
       95.45290,  82.73060,  71.83980,  62.52160,  54.56200,  47.90100,
       42.18700,  37.22830,  32.90820,  29.13450,  25.84930,  22.98280,
       20.47670,  18.24150,  16.24750,  14.49840,  12.95440,  11.55320,
       10.30320,   9.23414,   8.28228,   7.42266,   6.65575,   5.94844,
        5.40609,   4.88756,   4.45345,   4.03129,   3.68882,   3.34634};
  ///////////////////////////////

  // Assign absorption and refraction to materials

  // Quartz
  G4MaterialPropertiesTable* QuartzMPT = new G4MaterialPropertiesTable();
  QuartzMPT->AddProperty("RINDEX", PhotonEnergy, QuartzRefractiveIndex, num);
  QuartzMPT->AddProperty("ABSLENGTH", PhotonEnergy, QuartzAbsorption, num);

  // C4F10
  G4MaterialPropertiesTable* C4F10MPT = new G4MaterialPropertiesTable();
  C4F10MPT->AddProperty("RINDEX", PhotonEnergy, C4F10RefractiveIndex, num);
  C4F10MPT->AddProperty("ABSLENGTH", PhotonEnergy, Absorption, num);

  // Aerogel (LHCb-rich 1)
  G4MaterialPropertiesTable* AerogelMPT = new G4MaterialPropertiesTable();
  AerogelMPT->AddProperty("RINDEX", PhotonEnergy, AerogelRefractiveIndex, num);
  AerogelMPT->AddProperty("ABSLENGTH", PhotonEnergy, Absorption, num);
  
  // Cladding (of the fiber) only for the fiber aplication
  G4MaterialPropertiesTable* CladdingMPT = new G4MaterialPropertiesTable();
  CladdingMPT->AddProperty("RINDEX", PhotonEnergy, CladdingRefractiveIndex, num);
  CladdingMPT->AddProperty("ABSLENGTH", PhotonEnergy, QuartzAbsorption, num);

  // Air
  G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
  AirMPT->AddProperty("RINDEX", PhotonEnergy, AirRefractiveIndex, num);
  AirMPT->AddProperty("ABSLENGTH", PhotonEnergy, AirAbsorption, num);

  //Aerogel A Type A @ 420.0 nm --> 1.04611
  G4MaterialPropertiesTable* Aerogel_A_MPT = new G4MaterialPropertiesTable();
  Aerogel_A_MPT->AddProperty("RINDEX"  , PhotonEnergyAerogel, AerogelRefractiveIndexA, numAerogel);
  Aerogel_A_MPT->AddProperty("RAYLEIGH", PhotonEnergyAerogel, AerogelRayleighScattLA , numAerogel);

  //Aerogel B Type B @ 420.0 nm --> 1.05515
  G4MaterialPropertiesTable* Aerogel_B_MPT = new G4MaterialPropertiesTable();
  Aerogel_B_MPT->AddProperty("RINDEX"  , PhotonEnergyAerogel, AerogelRefractiveIndexB, numAerogel);
  Aerogel_B_MPT->AddProperty("RAYLEIGH", PhotonEnergyAerogel, AerogelRayleighScattLB , numAerogel);


  // Assign this material to the bars
  SiO2->SetMaterialPropertiesTable(QuartzMPT);
  SiO2_cladd->SetMaterialPropertiesTable(CladdingMPT);
  C4F10->SetMaterialPropertiesTable(C4F10MPT);
  Aerogel->SetMaterialPropertiesTable(AerogelMPT);
  AerogelA->SetMaterialPropertiesTable(Aerogel_A_MPT);
  AerogelB->SetMaterialPropertiesTable(Aerogel_B_MPT);
  Air->SetMaterialPropertiesTable(AirMPT);

}

void DetectorConstruction::ConstructField()
{
  //magnetic field
  static G4bool fieldIsInitialized = false;
  if(!fieldIsInitialized){
    G4FieldManager* fieldMgr
      = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
    fieldMgr->GetChordFinder()->SetDeltaChord(1.0*mm);
    fieldIsInitialized = true;    
  }
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{

  //xmlarichdataUtil();

  ConstructField();  

  G4double world_sizeX = 4.0*m;
  G4double world_sizeY = 4.0*m;
  G4double world_sizeZ = 6.0*m;

  // 
  // Define World Volume
  //
  G4VSolid *world_solid = new G4Box("World",world_sizeX/2.0,world_sizeY/2.0,world_sizeZ/2.0);
  G4LogicalVolume *world_logical = new G4LogicalVolume(world_solid,Air,"World");
  G4VPhysicalVolume *world_physical = new G4PVPlacement(0,G4ThreeVector(),world_logical,"World",0,false,0);

  G4double maxStep   = 0.1*mm;
  G4double maxLength = 2.0*m;
  G4double maxTime   = 20.0*ns; 
  G4double minEkin   = 1.0/100*MeV;
  G4double mionRang  = 0.01*mm;
  stepLimit = new G4UserLimits(maxStep,maxLength,maxTime,minEkin,mionRang);

  return world_physical;
}
