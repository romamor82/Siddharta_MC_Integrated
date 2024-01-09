#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaTrackerSD.h"
#include "../include/SiddhartaHisto.h"

#include <G4HCofThisEvent.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>

SiddhartaTrackerSD::SiddhartaTrackerSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname = "trackerCollection");
}

SiddhartaTrackerSD::~SiddhartaTrackerSD() {}

void SiddhartaTrackerSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new SiddhartaTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
  static G4int HCID = -1;

  if (HCID < 0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection(HCID, trackerCollection);
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    nHitFull = 0;
    for (G4int ichip=0; ichip<48; ichip++) {
      for (G4int isdd=0; isdd<8; isdd++) {
        sddNb[ichip*8 + isdd] = 0;
        chipNb[ichip*8 + isdd] = 0;
        nHit[ichip][isdd] = 0;
        sddEnergy[ichip][isdd] = 0.0;
      }
    }
    analysis->histo->ntuData.nHitSDD = 0;
    trackhit.clear();

    for (G4int i=0; i<200; i++) {
      for (G4int j=0; j<20; j++) {
        analysis->histo->ntuData.particleNameVertex[i][j] = '\0';
        analysis->histo->ntuData.materialVertex[i][j] = '\0';
      }
      analysis->histo->ntuData.matvert[i] = 0;
    }
  }
}

G4bool SiddhartaTrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    for (G4int i=0; i<200; i++) {
      analysis->histo->ntuData.matvert[i] = 0;
    }
    G4double edep = aStep->GetTotalEnergyDeposit();

    if (edep == 0.)
      return false;

    SiddhartaTrackerHit* newHit = new SiddhartaTrackerHit();
    newHit->SetTrackID(aStep->GetTrack()->GetTrackID());

    G4int copyNbSDD = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
    G4int copyNbSDDBOX = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(2)->GetCopyNo();
    G4int copyNbSDDCeramic = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1)->GetCopyNo();
    copyNbSDD = copyNbSDD + 4*(copyNbSDDCeramic - 1);

    newHit->SetSDDNb(copyNbSDD);
    newHit->SetCHIPNb(copyNbSDDBOX);
    newHit->SetEdep(edep);
    newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
    trackerCollection->insert(newHit);

    G4String pname = "";
    pname = aStep->GetTrack()->GetDefinition()->GetParticleName();
    G4String matname = "";
    matname = aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetMaterial()->GetName();
    sddEnergy[copyNbSDDBOX-1][copyNbSDD-1] += edep;

    if (nHit[copyNbSDDBOX-1][copyNbSDD-1] == 0) {
      nHit[copyNbSDDBOX-1][copyNbSDD-1]++;
      nHitFull++;
      sddNb[nHitFull-1] = copyNbSDD - 1;
      chipNb[nHitFull-1] = copyNbSDDBOX - 1;
      analysis->histo->ntuData.NoSDD[nHitFull-1] = 8*chipNb[nHitFull-1] + sddNb[nHitFull-1];
      trackhit[aStep->GetTrack()->GetTrackID()] = 8*(copyNbSDDBOX - 1) + copyNbSDD;

      analysis->histo->ntuData.XYZvertexSDD[nHitFull-1][0] = (aStep->GetTrack()->GetVertexPosition())[0]/mm;
      analysis->histo->ntuData.XYZvertexSDD[nHitFull-1][1] = (aStep->GetTrack()->GetVertexPosition())[1]/mm;
      analysis->histo->ntuData.XYZvertexSDD[nHitFull-1][2] = (aStep->GetTrack()->GetVertexPosition())[2]/mm;
      analysis->histo->ntuData.XYZSDD[nHitFull-1][0] = (aStep->GetTrack()->GetPosition())[0]/mm;
      analysis->histo->ntuData.XYZSDD[nHitFull-1][1] = (aStep->GetTrack()->GetPosition())[1]/mm;
      analysis->histo->ntuData.XYZSDD[nHitFull-1][2] = (aStep->GetTrack()->GetPosition())[2]/mm;

      G4int isdd = analysis->histo->ntuData.NoSDD[nHitFull-1];
      G4double pposx = (aStep->GetPreStepPoint()->GetPosition())[0]/mm;
      G4double sposx = analysis->SDD_pos[isdd][0]/mm;
      G4double pposy = (aStep->GetPreStepPoint()->GetPosition())[1]/mm;
      G4double sposy = analysis->SDD_pos[isdd][1]/mm;
      G4double pposz = (aStep->GetPreStepPoint()->GetPosition())[2]/mm;
      G4double sposz = analysis->SDD_pos[isdd][2]/mm;
      G4double th = -analysis->SDD_angle[isdd];
      G4double diffx =  cos(th)*(pposx - sposx) + sin(th)*(pposz - sposz);
      G4double diffy =  pposy - sposy;
      G4double diffz = -sin(th)*(pposx - sposx) + cos(th)*(pposz - sposz);
/*
      if (abs(diffx) > 0.5*450.*um) {
        G4cout << "Problem in X "<< G4endl;
        G4cout << "--- "<< isdd << G4endl;
        G4cout << "-X- "<< pposx << " " << sposx << " " << diffx << G4endl;
        G4cout << "-Y- "<< pposy << " " << sposy << " " << diffy << G4endl;
        G4cout << "-Z- "<< pposz << " " << sposz << " " << diffz << G4endl;
      }
      if (abs(diffz) > 0.5*1.*cm) {
        G4cout << "Problem in Z "<< G4endl;
        G4cout << "--- " << isdd << G4endl;
        G4cout << "-X- " << pposx << " "<< sposx << " " << diffx << G4endl;
        G4cout << "-Y- " << pposy << " "<< sposy << " " << diffy << G4endl;
        G4cout << "-Z- " << pposz << " "<< sposz << " " << diffz << G4endl;
      }
      if (abs(diffy) > 0.5*1.*cm) {
        G4cout << "Problem in  Y " << G4endl;
        G4cout << "--- " << isdd << G4endl;
        G4cout << "-X- " << pposx << " " << sposx << " " << diffx << G4endl;
        G4cout << "-Y- " << pposy << " " << sposy << " " << diffy << G4endl;
        G4cout << "-Z- " << pposz << " " << sposz << " " << diffz << G4endl;
      }
*/
      analysis->histo->ntuData.XYZInSDD[nHitFull-1][0] = diffx/mm;
      analysis->histo->ntuData.XYZInSDD[nHitFull-1][1] = diffy/mm;
      analysis->histo->ntuData.XYZInSDD[nHitFull-1][2] = diffz/mm;
      analysis->histo->ntuData.KinvertexSDD[nHitFull-1] = (aStep->GetTrack()->GetVertexKineticEnergy())/eV;
      analysis->histo->ntuData.parentID[nHitFull-1] = aStep->GetTrack()->GetParentID();
      analysis->histo->ntuData.TimeSDD[nHitFull-1] = aStep->GetTrack()->GetGlobalTime()/ns;
      analysis->histo->ntuData.particleVertexPDGEncoding[nHitFull-1] = aStep->GetTrack()->GetDefinition()->GetAntiPDGEncoding();
      int ilen = strlen(matname);

      for (G4int i=0; i<ilen; i++) {
        analysis->histo->ntuData.materialVertex[nHitFull-1][i] = (matname.data())[i];

        if ((i < 9) && ((matname.data())[i] < 123) && ((matname.data())[i] > 31)) {
          analysis->histo->ntuData.matvert[nHitFull-1]=analysis->histo->ntuData.matvert[nHitFull-1]*100 + ((matname.data())[i] - 30);
        }
      }
      for (G4int i=ilen; i<20; i++) {
        analysis->histo->ntuData.materialVertex[nHitFull-1][i] = '\0';
      }
      analysis->histo->ntuData.matVertIDFromMap[nHitFull-1] = analysis->histo->getMatIDFromString(matname);
    } else {
      if( trackhit[aStep->GetTrack()->GetTrackID()] == 8*(copyNbSDDBOX-1) + copyNbSDD ) {
      } else {
        if (aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetMaterial()->GetName() != "Silicon") {
          nHit[copyNbSDDBOX-1][copyNbSDD-1]++;
          trackhit[aStep->GetTrack()->GetTrackID()] = 8*(copyNbSDDBOX - 1) + copyNbSDD ;
        }
      }
    }
  }
  return true;
}

void SiddhartaTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    G4int ih = 0;
    analysis->histo->ntuData.nHitSDD = ih;

    for (G4int i=0; i<200; i++) {
      analysis->histo->ntuData.matvert[i] = 0;
      analysis->histo->ntuData.matVertIDFromMap[i] = 0;
    }

    if (nHitFull >= 700)
      G4cout << "nHitFull > 700 ? (" << nHitFull << ")" << G4endl;

    for (G4int index=0; index<nHitFull; index++) {
      if (analysis->YesHistos)
        analysis->histo->fillHisto("6", sddEnergy[chipNb[index]][sddNb[index]], 1.);
      ih++;
      analysis->histo->ntuData.nHitSDD = ih;
      analysis->histo->ntuData.NoSDD[ih-1] = analysis->histo->ntuData.NoSDD[index];
      analysis->histo->ntuData.XYZvertexSDD[ih-1][0] = analysis->histo->ntuData.XYZvertexSDD[index][0];
      analysis->histo->ntuData.XYZvertexSDD[ih-1][1] = analysis->histo->ntuData.XYZvertexSDD[index][1];
      analysis->histo->ntuData.XYZvertexSDD[ih-1][2] = analysis->histo->ntuData.XYZvertexSDD[index][2];
      analysis->histo->ntuData.KinvertexSDD[ih-1] = analysis->histo->ntuData.KinvertexSDD[index];
      analysis->histo->ntuData.parentID[ih-1] = analysis->histo->ntuData.parentID[index];
      analysis->histo->ntuData.TimeSDD[ih-1] = analysis->histo->ntuData.TimeSDD[index];
      analysis->histo->ntuData.particleVertexPDGEncoding[ih-1] = analysis->histo->ntuData.particleVertexPDGEncoding[index];
      int ilen = strlen(analysis->histo->ntuData.materialVertex[index]);

      for (G4int i=0; i<ilen+1; i++) {
        char chh = analysis->histo->ntuData.materialVertex[index][i];
        analysis->histo->ntuData.materialVertex[ih-1][i] = chh;

        if ((i < 11) && (chh > 31) && (chh < 123)) {
          analysis->histo->ntuData.matvert[ih-1]=analysis->histo->ntuData.matvert[ih-1]*100 + (chh - 30);
        }
      }
      analysis->histo->ntuData.matVertIDFromMap[ih-1] = analysis->histo->getMatIDFromString(analysis->histo->ntuData.materialVertex[index]);

      for(G4int i=ilen;i<20;i++) {
        analysis->histo->ntuData.materialVertex[ih-1][i] = '\0';
      }
      analysis->histo->ntuData.EnergySDD[ih-1] = sddEnergy[chipNb[index]][sddNb[index]]/eV;
      analysis->histo->ntuData.NbHitSDD[ih-1] = nHit[chipNb[index]][sddNb[index]];
    }
  }
}
