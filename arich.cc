#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4UI_USE_TCSH
#include "G4UItcsh.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "DetectorConstruction.hh"
#include "DetectorConstructionG1.hh"
#include "DetectorConstructionG2.hh"
#include "DetectorConstructionG3.hh"
#include "DetectorConstructionG4.hh"
#include "DetectorConstructionG5.hh"
#include "DetectorConstructionG6.hh"
#include "DetectorConstructionG7.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "StackingAction.hh"
#include "TrackingAction.hh"
//#include "fTOF_UImessenger.hh"
#include "SteppingVerbose.hh"

#include "Randomize.hh"

#include "G4SystemOfUnits.hh"

#include "TString.h"

//using namespace CLHEP;

int main(int argc, char** argv)
{

  // Choose the Random Engine
  //  HepRandom::setTheEngine(new RanecuEngine);
  // Seed the random number generator manually
  //

  if(argc!=12){
    G4cout<<" ERROR of the input parameters !!! "<<G4endl
	  <<"      [0] - vis.mac or run.mac or *.mac "<<G4endl
	  <<"      [1] - seed "<<G4endl
	  <<"      [2] - output root file name"<<G4endl
	  <<"      [3] - name of the particle ( e+, e-, mu+, mu-, pi+, pi-, kaon+, kaon-, proton, gamma, opticalphoton)"<<G4endl;
    G4cout<<"      [4] - particle x0 cm"<<G4endl
	  <<"      [5] - particle y0 cm"<<G4endl
	  <<"      [6] - particle z0 cm"<<G4endl;
    G4cout<<"      [7] - particle momentum (GeV/c) or optical photon wavelength (nm)"<<G4endl
	  <<"      [8] - particle theta (deg)"<<G4endl
	  <<"      [9] - particle phi (deg)"<<G4endl
	  <<"      [10] - geometry identification number"<<G4endl;
    return 0;    
  }
  else{
    G4cout<<"  Input parameters         "<<G4endl
	  <<"     mac file              "<<argv[1]<<G4endl
	  <<"     seed                  "<<argv[2]<<G4endl
	  <<"     output root file name "<<argv[3]<<G4endl
	  <<"     name of the particle  "<<argv[4]<<G4endl;
    G4cout<<"     particle x0           "<<argv[5]<<" cm"<<G4endl
	  <<"     particle y0           "<<argv[6]<<" cm"<<G4endl
	  <<"     particle z0           "<<argv[7]<<" cm"<<G4endl;
    G4cout<<"     particle momentum     "<<argv[8]<<" GeV/c"<<G4endl
	  <<"     particle theta        "<<argv[9]<<" deg"<<G4endl
	  <<"     particle phi          "<<argv[10]<<" deg"<<G4endl
	  <<"     geometry ID           "<<argv[11]<<G4endl;
  }

  G4long myseed = 345354;
  myseed = atoi(argv[2]);
  G4cout<<" myseed - "<<myseed<<G4endl;

  CLHEP::HepRandom::setTheSeed(myseed);

  // Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // Construct mandatory initialization classes
  G4VUserDetectorConstruction* detector;
  int geomID = atoi(argv[11]);
  if(geomID == 0){
    detector = new DetectorConstruction;
  }
  else if(geomID == 1){
    detector = new DetectorConstructionG1;
  }
  else if(geomID == 2){
    detector = new DetectorConstructionG2;
  }
  else if(geomID == 3){
    detector = new DetectorConstructionG3;
  }
  else if(geomID == 4){
    detector = new DetectorConstructionG4;
  }
  else if(geomID == 5){
    detector = new DetectorConstructionG5;
  }
  else if(geomID == 6){
    detector = new DetectorConstructionG6;
  }
  else if(geomID == 7){
    detector = new DetectorConstructionG7;
  }
  else{
    assert(0);
  }
  runManager->SetUserInitialization(detector);
  G4VUserPhysicsList* physics = new PhysicsList;
  runManager->SetUserInitialization(physics);

  // Construct User Action Classes
  RunAction* runAction = new RunAction;
  runManager->SetUserAction(runAction);
  //
  G4String rootFileName = argv[3];
  G4cout<<rootFileName<<G4endl;
  runAction->SetOutputFileName(rootFileName);

  PrimaryGeneratorAction* genAction = 
    new PrimaryGeneratorAction();
  runManager->SetUserAction(genAction);
  G4String particleName = argv[4];
  G4double x0 = atof(argv[5]);
  G4double y0 = atof(argv[6]);
  G4double z0 = atof(argv[7]);
  G4double particleMomentum = atof(argv[8]);
  G4double particleTheta = atof(argv[9]);
  G4double particlePhi = atof(argv[10]);

  G4cout<<"particle Name = "<<particleName<<G4endl;
  G4cout<<"x0 = "<<x0<<" cm"<<G4endl;
  G4cout<<"y0 = "<<y0<<" cm"<<G4endl;
  G4cout<<"z0 = "<<z0<<" cm"<<G4endl;
  G4cout<<"particle Momentum = "<<particleMomentum<<G4endl;
  G4cout<<"particle Theta = "<<particleTheta<<G4endl;
  G4cout<<"particle Phi = "<<particlePhi<<G4endl;
  genAction->SetParticleName(particleName);
  genAction->SetXYZvertex( x0, y0, z0);
  if(particleName == "opticalphoton")
    genAction->SetParticleMomentum(particleMomentum);
  else
    genAction->SetParticleMomentum(particleMomentum*GeV);
  genAction->SetThetaAngle(particleTheta*deg);
  genAction->SetPhiAngle(particlePhi*deg);

  TString ffName = "SteppingAction_";
  TString fffName = rootFileName;
  ffName += fffName;
  ffName += ".txt";
  G4cout<<"ffName = "<<ffName<<G4endl;
  SteppingAction* steppingAction = new SteppingAction(genAction,ffName);
  runManager->SetUserAction(steppingAction);

  EventAction* eventAction = new EventAction(runAction,
							 steppingAction);
  runManager->SetUserAction(eventAction);
  StackingAction *stackingAction = new StackingAction;
  runManager->SetUserAction(stackingAction);
  runManager->SetUserAction(new TrackingAction);

  eventAction->SetStackingAction(stackingAction);
  eventAction->SetPrimGenerator(genAction);
  // Setup to be able to define some custom commands;
  //fTOF_UImessenger* messenger = new fTOF_UImessenger(runAction, genAction);

  // Initialize G4 kernel
  runManager->Initialize();

  // Get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (argc == 1) {  // Define UI terminal for interactive mode
    UI->ApplyCommand("/run/verbose 1");
    UI->ApplyCommand("/event/verbose 1");
    UI->ApplyCommand("/tracking/verbose 1");

    //////////
    //G4VisManager* visManager = new G4VisExecutive;
    //visManager->Initialize();
    /////////

#if defined(G4UI_USE_TCSH)
    G4UIsession* session = new G4UIterminal(new G4UItcsh);
#else
    G4UIsession* session = new G4UIterminal;
#endif
    session->SessionStart();

    delete session;
  }
  else {

    G4String fileName = argv[1];
    if(fileName.contains("vis")){
      G4cout<<"VIS fileName "<<fileName<<G4endl;
#ifdef G4VIS_USE
      // Visualization manager
      G4VisManager* visManager = new G4VisExecutive;
      visManager->Initialize();
#endif    
      
      G4String command = "/control/execute ";
      UI->ApplyCommand(command+fileName);
      
#ifdef G4VIS_USE
      delete visManager;
#endif
    }
    else{
      G4cout<<"NO VIS fileName "<<fileName<<G4endl;
      G4String command = "/control/execute ";
      UI->ApplyCommand(command+fileName);    
    }

  }

  // Start a run
  //  G4int numberOfEvent = 3;
  //  runManager->BeamOn(numberOfEvent);

  // Job termination
  //
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  //delete messenger;
  //delete visManager;
  //G4cout<<"  delete runManager "<<G4endl;
  delete runManager;
  
  return 0;
}
