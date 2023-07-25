#include "../include/SiddhartaPrimaryGeneratorAction.h"
#include "../include/SiddhartaDetectorConstruction.h"
#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaEventAction.h"
#include "../include/SiddhartaRunAction.h"
#include "../include/PhiMassFunction.h"
#include "../include/SiddhartaHisto.h"
#include "../include/SiddhartaCard.h"

#include <G4VShortLivedParticle.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4ParticleTable.hh>
#include <G4LorentzVector.hh>
#include <G4ExcitedMesons.hh>
#include <G4ParticleGun.hh>
#include <G4RunManager.hh>
#include <Randomize.hh>
#include <globals.hh>
#include <G4Event.hh>

#include <TF1.h>

#include <iostream>
#include <fstream>

int Sign (double a) 
{
  if (a > 0)
    return +1;
  else if (a < 0)
    return -1;
  else
    return 0;
}

SiddhartaPrimaryGeneratorAction::SiddhartaPrimaryGeneratorAction(SiddhartaDetectorConstruction* myDC) : myDetector(myDC)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleEnergy(0.0);
  SiddhartaCard* mycard = SiddhartaCard::getInstance();

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4int theBeam = mycard->variables["Beam"];

  if (theBeam == 0) {
    particleVertex = particleTable->FindParticle("e-");
  } else if (theBeam == 1) {
    particleVertex = particleTable->FindParticle("phi");
  } else if (theBeam == 2) {
    particleVertex = particleTable->FindParticle("kaon-");
  } else if (theBeam == 3) {
    particleVertex = particleTable->FindParticle("gamma");
  }

  if (particleVertex)
    particleGun->SetParticleDefinition(particleVertex);
}

SiddhartaPrimaryGeneratorAction::~SiddhartaPrimaryGeneratorAction()
{
  delete particleGun;
}

void SiddhartaPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  SiddhartaCard* mycard = SiddhartaCard::getInstance();
  G4ThreeVector vertex;
  G4double Kinener;
  SiddhartaAnalysisManager* analysis;

  if (particleVertex->GetParticleName() == "phi") {
    TF1* phiMass = new TF1("phiMass", phiMassFunc, 1000, 1100., 2);
    phiMass->FixParameter(0, mphi_pdg);
    phiMass->FixParameter(1, phi_Gm);
    phiMass->SetNpx(10000);
    G4double eMass = 0.511*MeV;
    G4double eBeam = 510.0*MeV;
    G4double eSigma = 0.2*MeV;
    G4double halfCrossingAngle = 25.*mrad;
    G4double eElectron = 0;
    G4double ePositron = 0;
    G4double pElectron = 0;
    G4double pPositron = 0;
    G4LorentzVector pEl;
    G4LorentzVector pPo;
    G4LorentzVector pPhi;
    G4double pMom = 0;
    G4double invMphi = 0;
    G4double prob = 0;
    G4int loopCounter = 0;

    while (1) {
      G4double eSmear = G4RandGauss::shoot(0., 1.);
      G4double pSmear = G4RandGauss::shoot(0., 1.);

      eElectron = eBeam + eSmear*eSigma;
      ePositron = eBeam + pSmear*eSigma;
      pElectron = sqrt(eElectron*eElectron - eMass*eMass);
      pPositron = sqrt(ePositron*ePositron - eMass*eMass);

      pEl[3] = eElectron;
      pEl[0] = pElectron*sin(halfCrossingAngle);
      pEl[1] = 0.0;
      pEl[2] = pElectron*cos(halfCrossingAngle);

      pPo[3] = ePositron;
      pPo[0] = pPositron*sin(halfCrossingAngle);
      pPo[1] = 0.0;
      pPo[2] = -pPositron*cos(halfCrossingAngle);

      pPhi[0] = pEl[0] + pPo[0];
      pPhi[1] = pEl[1] + pPo[1];
      pPhi[2] = pEl[2] + pPo[2];
      pPhi[3] = pEl[3] + pPo[3];

      pMom = sqrt(pPhi[0]*pPhi[0]+pPhi[1]*pPhi[1]+pPhi[2]*pPhi[2]);
      invMphi = pPhi.m();
      prob = phiMass->Eval(invMphi)/phiMass->Eval(mphi_pdg);

      if (G4UniformRand() < prob) {
        break;
      } else if (loopCounter == 1000) {  // Set some variable here {
        G4cout << "----------------------------------------------------------"
               << "------ Looped for 1000 times without phi production ------"
               << "----------------------------------------------------------" << G4endl;
        break;
      }
      loopCounter ++ ;
    }
    particleGun->SetParticleTime(0.0*ns);
    particleGun->SetParticleMomentum(pMom);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(pPhi[0],pPhi[1],pPhi[2]));
    particleGun->SetParticlePolarization(G4ThreeVector(0.,0.,1.));
    G4ThreeVector bSigma = G4ThreeVector(2*mm,20*um,1.0*cm);

    for (G4int i = 0; i<3; i++) {
      vertex[i] = G4RandGauss::shoot(0.0, bSigma[i]);
    }
    particleGun->SetParticlePosition(vertex);
  } else if (particleVertex->GetParticleName() == "e-" ) {
    particleGun->SetParticleTime(0.0*ns);
    G4double ranflat = G4UniformRand();
    G4double ranflat1 = G4UniformRand();

    float nbOfSeconds = mycard->variables["nbOfSeconds"];
    float beamCurrent = mycard->variables["beamCurrent"];

    if (mycard->GetNoDataStatus() == true) {
      if (mycard->variables["externalBeam"] == true) {
        G4cout << "Error opening beam file, using internal generator!" << G4endl;
      }
      G4int updown = (ranflat1>.5)*2 - 1;
      vertex[0] = 0.0*mm;
      vertex[1] = updown*27.5*mm;
      vertex[2] = (((ranflat - .5)*240.) + 470. + 120)*mm;
      particleGun->SetParticleMomentum(510.0*MeV);
      particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.1*updown,-.99));
      particleGun->SetParticlePosition(vertex);
    } else {
      G4double ang_out = 4.868*deg ;
      G4double xa = 17.7 * mm ;
      G4double za = 535. * mm ;
      G4double ang_in = atan(xa/za);
      G4double AO = sqrt(xa*xa + za*za);

      if (mycard->GetInteractive() == false) {
        if (mycard->GetDoneStatus() == 0) {
          do {
            if (!mycard->BeamEof()) {
			  int cnt = mycard->GetBeam();
			  xx=mycard->xx;
              xxp=mycard->xxp;
              yy=mycard->yy;
              yyp=mycard->yyp;
              zz=mycard->zz;
              de=mycard->de;
              rate=mycard->rate;
              turn=mycard->turn;
	      if ((cnt%1000)==0) {
                G4cout << " beam sample. db.recno=" << cnt << " par: " << xx << " " << xxp << " " << yy << " " << yyp << " " << zz << " " << de << " " << rate << " " << turn << "--- \n";
              }
              mycard->SetEventMult(rate*0.1*nbOfSeconds*beamCurrent); //2023 Feb 16 this was wrongly defined here, has to be remembered for eventmul calls
            } else {
              G4cout << "Beam data abnormally ended, exiting! \n";
              G4RunManager::GetRunManager()->AbortRun();  	//2023 Feb 16 exit does not work, err mess are filling the disk- investigate
            }
          } while (mycard->GetEventMult() < 1);
        }
        if (abs(zz*m) > AO) {
          vertex[0] = xa*Sign(zz) + (zz*m - AO*Sign(zz))*sin(ang_out) + xx*m*cos(ang_out);
          vertex[1] = yy*m + (ranflat1 - .5)*mm;
          vertex[2] = za*Sign(zz) + (zz*m - AO*Sign(zz))*cos(ang_out) - xx*m*sin(ang_out) + (ranflat - .5)*mm;
        } else {
          vertex[0] = zz*m*sin(ang_in) + xx*m*cos(ang_in);
          vertex[1] = yy*m + (ranflat1 - .5)*mm;
          vertex[2] = zz*m*cos(ang_in) - xx*m*sin(ang_in) + (ranflat - .5)*mm;
        }
        particleGun->SetParticleMomentum((1 + de)*510.0*MeV);
        particleGun->SetParticleMomentumDirection(G4ThreeVector(xxp + ang_out, yyp,.99));
        particleGun->SetParticlePosition(vertex);
        mycard->AddDoneByOne();

        if (mycard->GetDoneStatus() == mycard->GetEventMult()) {
          mycard->SetDoneStatus(0);
        }
      } else {
        if (!mycard->BeamEof()) {
          mycard->GetBeam();
		  xx=mycard->xx;
          xxp=mycard->xxp;
          yy=mycard->yy;
          yyp=mycard->yyp;
          zz=mycard->zz;
          de=mycard->de;
          rate=mycard->rate;
          turn=mycard->turn;
          if (abs(zz*m) > AO) {
            vertex[0] = xa*Sign(zz) + (zz*m - AO*Sign(zz))*sin(ang_out) + xx*m*cos(ang_out);
            vertex[1] = yy*m + (ranflat1 - .5)*mm;
            G4cout << " out_case:  \n";
            vertex[2] = za*Sign(zz) + (zz*m - AO*Sign(zz))*cos(ang_out) - xx*m*sin(ang_out) + (ranflat - .5)*mm;
          } else {
            vertex[0] = zz*m*sin(ang_in) + xx*m*cos(ang_in);
            vertex[1] = yy * m + (ranflat1 - .5)*mm;
            vertex[2] = zz*m*cos(ang_in) - xx*m*sin(ang_in) + (ranflat - .5)*mm;
          }
          particleGun->SetParticleMomentum((1 + de)*510.0*MeV);
          particleGun->SetParticleMomentumDirection(G4ThreeVector(xxp + ang_out, yyp,.99));
          particleGun->SetParticlePosition(vertex);
          G4cout << "data:    " << xx << " " << xxp << " " << yy << " " << yyp << " " << zz << " " << de << " " << rate << " " << turn << "\n";
          G4cout << "ang_in= " << ang_in << "  AO=" << AO << "  x=" << vertex[0] << " z=" << vertex[2] << "\n";
        } else {
//          nodata = true;
          mycard->SetNoDataStatus(true);
          G4cout << "Beam file data exausted, using internal generator! \n";
          G4int updown = (ranflat1 > .5)*2 - 1;
          vertex[0] = 0.0*mm ;
          vertex[1] = updown*27.5*mm ; //signof(ranflat1-.5)...
          vertex[2] = (((ranflat - .5) * 240.) + 470. + 120)* mm ;
          particleGun->SetParticleMomentum(510.0*MeV);
          particleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0.1*updown, -.99));
          particleGun->SetParticlePosition(vertex);
        }
      }
    }
  } else if (particleVertex->GetParticleName() == "kaon-") {
    particleGun->SetParticleTime(0.0*ns);
    vertex[0] = 0;
    vertex[1] = 0.*cm;
    vertex[2] = 0.*cm;
    particleGun->SetParticleMomentum(126.35*MeV);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,1.,0.));
    particleGun->SetParticlePosition(vertex);
  } else if (particleVertex->GetParticleName() == "gamma") {
    particleGun->SetParticleTime(0.0*ns);
    vertex[0] = 0;
    vertex[1] = 0.*cm;
    vertex[2] = 0.*cm;
    particleGun->SetParticleEnergy(G4RandGauss::shoot(25., 5.) * keV);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,1.,0.));
    particleGun->SetParticlePosition(vertex);
  }
  analysis = SiddhartaAnalysisManager::getInstance();
  Kinener = particleGun->GetParticleEnergy();

  if (analysis->YesHistos) {
    analysis->histo->fillHisto("0", vertex[0], 1.);
    analysis->histo->fillHisto("1", vertex[1], 1.);
    analysis->histo->fillHisto("2", vertex[2], 1.);
    analysis->histo->fillHisto("3", Kinener, 1.);
    analysis->histo->fillHisto3("4", vertex[2], vertex[0], vertex[1], 1.);
    analysis->histo->ntuData.XYZvertex[0] = vertex[0];
    analysis->histo->ntuData.XYZvertex[1] = vertex[1];
    analysis->histo->ntuData.XYZvertex[2] = vertex[2];
    analysis->histo->ntuData.VertexMomentum[0] = particleGun->GetParticleMomentum()*(particleGun->GetParticleMomentumDirection())[0]/MeV;
    analysis->histo->ntuData.VertexMomentum[1] = particleGun->GetParticleMomentum()*(particleGun->GetParticleMomentumDirection())[1]/MeV;
    analysis->histo->ntuData.VertexMomentum[2] = particleGun->GetParticleMomentum()*(particleGun->GetParticleMomentumDirection())[2]/MeV;
  }
  particleGun->GeneratePrimaryVertex(anEvent);
}
