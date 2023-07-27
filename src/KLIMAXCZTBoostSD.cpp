#include "../include/SiddhartaAnalysisManager.h"
#include "../include/KLIMAXCZTBoostSD.h"
#include "../include/SiddhartaHisto.h"

#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>


KLIMAXCZTBoostSD::KLIMAXCZTBoostSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollectionKLCZTBoost");
}

KLIMAXCZTBoostSD::~KLIMAXCZTBoostSD() {}

void KLIMAXCZTBoostSD::Initialize(G4HCofThisEvent* HCE)
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
    TimeKLCZTBoost = -9999.;
    X = -99999.;
    Y = -99999.;
    Z = -99999.;
    analysis->histo->ntuData.TotalEnergyDepKLCZTBoost = -1000000. / eV;
    analysis->histo->ntuData.nHitCZTBoost = 0;

    for (int i=0; i<MaxCZTHits; i++) {
      analysis->histo->ntuData.pdgcodeKLCZTBoost[i] = -1000000;
      analysis->histo->ntuData.EnergyDepKLCZTBoost[i] = -1000000. / eV;
      analysis->histo->ntuData.TimeKLCZTBoost[i] = -1000000. / ns;
      analysis->histo->ntuData.XYZKLCZTBoost[0][i] = -1000000. / mm;
      analysis->histo->ntuData.XYZKLCZTBoost[1][i] = -1000000. / mm;
      analysis->histo->ntuData.XYZKLCZTBoost[2][i] = -1000000. / mm;
      analysis->histo->ntuData.kaonKinEKLCZTBoost[i] = -1000000. / eV;
    }
    analysis->histo->ntuData.gammaKinEKLCZTBoost = -1000000. / eV;
    gammacounter = 0;
  }
}

G4bool KLIMAXCZTBoostSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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
      G4cout << gammacounter << " photon in CZT Boost: " << aStep->GetTrack()->GetKineticEnergy() / eV  << G4endl;

      if (gammacounter == 1)
        analysis->histo->ntuData.gammaKinEKLCZTBoost = aStep->GetTrack()->GetKineticEnergy() / eV;
    }

    if (edep > sciEnergyMax) {
      sciEnergyMax = edep;
    }
    analysis->histo->ntuData.TotalEnergyDepKLCZTBoost = sciEnergy / eV;
    G4int CZTBoostcounter = analysis->histo->ntuData.nHitCZTBoost;
    TimeKLCZTBoost = aStep->GetTrack()->GetGlobalTime();
    X = (aStep->GetTrack()->GetPosition())[0];
    Y = (aStep->GetTrack()->GetPosition())[1];
    Z = (aStep->GetTrack()->GetPosition())[2];

    analysis->histo->ntuData.pdgcodeKLCZTBoost[CZTBoostcounter] = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();
    analysis->histo->ntuData.EnergyDepKLCZTBoost[CZTBoostcounter] = edep / eV;
    analysis->histo->ntuData.TimeKLCZTBoost[CZTBoostcounter] = TimeKLCZTBoost / ns;
    analysis->histo->ntuData.XYZKLCZTBoost[0][CZTBoostcounter] = X / mm;
    analysis->histo->ntuData.XYZKLCZTBoost[1][CZTBoostcounter] = Y / mm;
    analysis->histo->ntuData.XYZKLCZTBoost[2][CZTBoostcounter] = Z / mm;
    analysis->histo->ntuData.kaonKinEKLCZTBoost[CZTBoostcounter] = aStep->GetTrack()->GetKineticEnergy() / eV;
    analysis->histo->ntuData.nHitCZTBoost = CZTBoostcounter + 1;
  }

  return true;
}

void KLIMAXCZTBoostSD::EndOfEvent(G4HCofThisEvent*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  gammacounter = 0;
  kaonCounter = 0;
  kaonKinE = 0;
}
