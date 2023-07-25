#include "include/SiddhartaPrimaryGeneratorAction.h"
#include "include/SiddhartaDetectorConstruction.h"
#include "include/SiddhartaTrackingAction.h"
#include "include/SiddhartaSteppingAction.h"
#include "include/SiddhartaPhysicsList.h"
#include "include/SiddhartaEventAction.h"
#include "include/SiddhartaRunAction.h"
#include "include/SiddhartaCard.h"

#include <G4TrajectoryDrawByParticleID.hh>
#include <G4TrajectoryDrawByCharge.hh>
#include <G4SystemOfUnits.hh>
#include <G4RunManager.hh>
#include <G4UImanager.hh>

#include <iostream>
#include <fstream>
#include <map>

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "SiddhartaSteppingVerbose.h"


int main(int argc, char** argv)
{

/*
Structure for external arguments - detector setup, type of the beam... -> To be modified with messenger objects!!
Idea for main:
- Moving the reading of the simulation parameters by messenger -> Executed by action files
- Executing;
 + RunManager
 + Constuction
 + Actions
 + Physics
 + UI
 + Seed (Initialization + saving?)
- Removing drawing and draw options to other class/structure
  */

  SiddhartaCard *mycard = SiddhartaCard::getInstance();
  if (argc > 2) {
    std::string cardFileName = argv[2];
    mycard->SetCardFileName(cardFileName);
  }
  mycard->Initialization();
  mycard->ReadCard();

/*
Setting different output level. Is it neccessary? Maybe moving to lower level.
 */
  G4VSteppingVerbose* verbosity = new SiddhartaSteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);

  G4RunManager *runManager = new G4RunManager;
  SiddhartaDetectorConstruction* detector = new SiddhartaDetectorConstruction;

  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new SiddhartaPhysicsList);

  runManager->SetUserAction(new SiddhartaPrimaryGeneratorAction(detector));
  runManager->SetUserAction(new SiddhartaRunAction);
  runManager->SetUserAction(new SiddhartaEventAction);
  runManager->SetUserAction(new SiddhartaTrackingAction);
  runManager->SetUserAction(new SiddhartaSteppingAction);

  runManager->Initialize();

  G4VisManager *visManager = new G4VisExecutive;
  visManager->Initialize();

  G4TrajectoryDrawByParticleID* drawManger = new G4TrajectoryDrawByParticleID;

  drawManger->SetDefault("cyan");
  drawManger->Set("kaon-","red");
  drawManger->Set("kaon+","red");
  drawManger->Set("mu-","yellow");
  drawManger->Set("mu+","yellow");
  drawManger->Set("pi-","green");
  drawManger->Set("pi+","green");
  drawManger->Set("pi0","white");
  drawManger->Set("anti_nu_mu",G4Colour(0.5, 0.5, 0.5));
  drawManger->Set("nu_mu",G4Colour(0.5, 0.5, 0.5));
  drawManger->Set("anti_nu_e",G4Colour(0.5, 0.5, 0.5));
  drawManger->Set("nu_e",G4Colour(0.5, 0.5, 0.5));
  drawManger->Set("e-","magenta");
  drawManger->Set("e+","magenta");
  drawManger->Set("gamma",G4Colour(0.5, 0.2, 0.8));

  visManager->RegisterModel(drawManger);

  G4cout << "visual present" << G4endl;

/*
?? Is it neccessary? Just maybe passing empty line or more details in the output  2023 Feb 11 no,it was a video driver test
 */
  //G4cout << "after===================================" << G4endl;

  G4UImanager *UImanager = G4UImanager::GetUIpointer();

  bool externalBeam = mycard->variables["externalBeam"];
  G4int theBeam = mycard->variables["Beam"];
  int howmany = 0;
  mycard->SetNoDataStatus(true);
  
  std::ifstream inputbeam;
  if (theBeam == 0 && externalBeam == true) {
    G4int backgroundType = mycard->variables["backgroundType"];
    char *infile = "touschek.dat";

/* No other background source file -> Possible simplification of these part of the code by separate class/structure, so the
code in the main is cleaner!! -> Losing up the tables -> vectors with safe use */ 

//2023 Feb 11 No, there are actually 3 background files, as in the orig code, to be dealt with separately, according to DAFNE vacuum conditions
    
    G4cout << "Opening background file : "<< infile << G4endl;
    mycard->OpenBeam(infile);

    float nbOfSeconds = mycard->variables["nbOfSeconds"];
    float beamCurrent = mycard->variables["beamCurrent"];
    int records = 0;

    if(!mycard->BeamOpen()) {
      G4cout << "Error opening beam file" << G4endl;
    } else {
      G4cout << "Opening beam file" << G4endl;
      mycard->SetNoDataStatus(false);
      double xx, xxp, yy, yyp, zz, de, rate, turn;
      while (!mycard->BeamEof()) {
		mycard->GetBeam();
        records++;
		xx=mycard->xx;
        xxp=mycard->xxp;
        yy=mycard->yy;
        yyp=mycard->yyp;
        zz=mycard->zz;
        de=mycard->de;
        rate=mycard->rate;
        turn=mycard->turn;
        howmany += rate*0.1*nbOfSeconds*beamCurrent; // the file provides information for 1sec with current = 10mA -> Why 0.1? Because of the 10 mA.
		if (records%10000==0) {
          G4cout << " " << xx << " " << xxp << " " << yy << " " << yyp << " " << zz << " " << de << " " << rate << " " << turn << " \n";
        }
      }
	  mycard->CloseBeam();
	  mycard->OpenBeam(infile);
	  mycard->ResetCounter();
	  mycard->SetDoneStatus(0);
	  mycard->SetEventMult(0);
 //     inputbeam.open(infile); //2023 Feb 11 this should not be comented, is a rewind and the only open during run
      G4cout << ":::::::::::: " << howmany << " events are to be generated for " << nbOfSeconds << " seconds of "
             << beamCurrent << " mA run ::::::::::::::::::::::::::::::::::::::::\n";
    }
  }

  if (argc != 1) {
    int numberOfEvents = atoi(argv[1]);
    if (numberOfEvents > 0) {
      if (theBeam == 0 && mycard->GetNoDataStatus() == false) {
//      done = 0;
        mycard->SetDoneStatus(0);
        runManager->BeamOn(howmany-1);
      } else {
        runManager->BeamOn(numberOfEvents);
      }
    } else {
      mycard->SetInteractiveTrue();
      G4UIExecutive *ui = new G4UIExecutive(argc, argv);
      UImanager->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();

      delete ui;
      delete visManager;
    }
  } else if (argc == 1) {
    mycard->SetInteractiveTrue();

    G4UIExecutive *ui = new G4UIExecutive(argc, argv);

    UImanager->ApplyCommand("/control/execute vis.mac");

    ui->SessionStart();

    delete ui;
    delete visManager;

  }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  if (theBeam == 0 && mycard->GetNoDataStatus() == false) {
    inputbeam.close();
  }

  delete runManager;
  return 0;
}
