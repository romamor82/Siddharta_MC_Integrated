#include "../include/SiddhartaAnalysisManager.h"
#include "../include/KLIMAXHPGeBoostSD.h"
#include "../include/SiddhartaHisto.h"

#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>


KLIMAXHPGeBoostSD::KLIMAXHPGeBoostSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollectionKLHPGeBoost");
}

KLIMAXHPGeBoostSD::~KLIMAXHPGeBoostSD() {}

void KLIMAXHPGeBoostSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new SiddhartaTrackerHitsCollection(SensitiveDetectorName,collectionName[0]);
  static G4int HCID = -1;

  if (HCID<0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection(HCID, trackerCollection);

  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    sciEnergy = 0.0;
    sciEnergyMax = 0.0;
    TimeKLHPGeBoost = -9999.;
    X = -99999.;
    Y = -99999.;
    Z = -99999.;
    analysis->histo->ntuData.TotalEnergyDepKLHPGeBoost = -1000000. / eV;
    analysis->histo->ntuData.nHitHPGeBoost = 0;

    for (int i=0; i<MaxCZTHits; i++) {
      analysis->histo->ntuData.pdgcodeKLHPGeBoost[i] = -1000000;
      analysis->histo->ntuData.EnergyDepKLHPGeBoost[i] = -1000000. / eV;
      analysis->histo->ntuData.TimeKLHPGeBoost[i] = -1000000. / ns;
      analysis->histo->ntuData.XYZKLHPGeBoost[0][i] = -1000000. / mm;
      analysis->histo->ntuData.XYZKLHPGeBoost[1][i] = -1000000. / mm;
      analysis->histo->ntuData.XYZKLHPGeBoost[2][i] = -1000000. / mm;
      analysis->histo->ntuData.kaonKinEKLHPGeBoost[i] = -1000000. / eV;
    }
    analysis->histo->ntuData.gammaKinEKLHPGeBoost = -1000000. / eV;
    gammacounter = 0;
  }
}

G4bool KLIMAXHPGeBoostSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep == 0.)
    return false;

  SiddhartaTrackerHit* newHit = new SiddhartaTrackerHit();
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  trackerCollection->insert(newHit);

  G4String pname = "";
  pname = aStep->GetTrack()->GetDefinition()->GetParticleName();
  G4String prematname = "";
  G4String prevolname = "";
  G4String postmatname = "";
  G4String postvolname = "";

  G4StepStatus prestepState = aStep->GetPreStepPoint()->GetStepStatus();
  G4StepStatus poststepState = aStep->GetPostStepPoint()->GetStepStatus();
  prematname = aStep->GetPreStepPoint()->GetMaterial()->GetName();
  prevolname = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  postmatname = aStep->GetPostStepPoint()->GetMaterial()->GetName();
  postvolname = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
  G4double Kinener = aStep->GetTrack()->GetMomentum().mag();
  G4ThreeVector pos = aStep->GetTrack()->GetPosition();
  G4String preproc = "undefined";
  G4String postproc = "undefined";

  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  sciEnergy += edep;

  if (analysis->YesHistos) {
    if (aStep->GetTrack()->GetDynamicParticle()->GetPDGcode() == 22) {
      gammacounter++;
//      G4cout << gammacounter << " photon in CZT Boost: " << aStep->GetTrack()->GetKineticEnergy() / eV  << G4endl;

      if (gammacounter == 1) {analysis->histo->ntuData.gammaKinEKLHPGeBoost = aStep->GetTrack()->GetKineticEnergy() / eV;}
    }

    if (edep > sciEnergyMax) {
      sciEnergyMax = edep;
    }
    analysis->histo->ntuData.TotalEnergyDepKLHPGeBoost = sciEnergy / eV;
    G4int HPGeBoostcounter = analysis->histo->ntuData.nHitHPGeBoost;
    TimeKLHPGeBoost = aStep->GetTrack()->GetGlobalTime();
    X = (aStep->GetTrack()->GetPosition())[0];
    Y = (aStep->GetTrack()->GetPosition())[1];
    Z = (aStep->GetTrack()->GetPosition())[2];

    analysis->histo->ntuData.pdgcodeKLHPGeBoost[HPGeBoostcounter] = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();
    analysis->histo->ntuData.EnergyDepKLHPGeBoost[HPGeBoostcounter] = edep / eV;
    analysis->histo->ntuData.TimeKLHPGeBoost[HPGeBoostcounter] = TimeKLHPGeBoost / ns;
    analysis->histo->ntuData.XYZKLHPGeBoost[0][HPGeBoostcounter] = X / mm;
    analysis->histo->ntuData.XYZKLHPGeBoost[1][HPGeBoostcounter] = Y / mm;
    analysis->histo->ntuData.XYZKLHPGeBoost[2][HPGeBoostcounter] = Z / mm;
    analysis->histo->ntuData.kaonKinEKLHPGeBoost[HPGeBoostcounter] = aStep->GetTrack()->GetKineticEnergy() / eV;
    analysis->histo->ntuData.nHitHPGeBoost = HPGeBoostcounter + 1;
  }

  return true;
}

void KLIMAXHPGeBoostSD::EndOfEvent(G4HCofThisEvent*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  gammacounter = 0;
  kaonCounter = 0;
  kaonKinE = 0;
}
