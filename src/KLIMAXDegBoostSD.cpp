#include "../include/SiddhartaAnalysisManager.h"
#include "../include/KLIMAXDegBoostSD.h"
#include "../include/SiddhartaHisto.h"

#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>


KLIMAXDegBoostSD::KLIMAXDegBoostSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollectionKLDBoost");
}

KLIMAXDegBoostSD::~KLIMAXDegBoostSD() {}

void KLIMAXDegBoostSD::Initialize(G4HCofThisEvent* HCE)
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
    TimeKLDBoost = -9999.;
    X = -99999.;
    Y = -99999.;
    Z = -99999.;
    analysis->histo->ntuData.KLDBooststop[0] = -1000000 / mm;
    analysis->histo->ntuData.KLDBooststop[1] = -1000000 /  mm;
    analysis->histo->ntuData.KLDBooststop[2] = -1000000 / mm;
    analysis->histo->ntuData.EnergyDepKLDBoost = -1000000 / eV;
    analysis->histo->ntuData.TimeKLDBoost = -1000000 / ns;
    analysis->histo->ntuData.XYZKLDBoost[0] = -1000000 / mm;
    analysis->histo->ntuData.XYZKLDBoost[1] = -1000000 / mm;
    analysis->histo->ntuData.XYZKLDBoost[2] = -1000000 / mm;
    analysis->histo->ntuData.kaonKinEKLDBoost = -1000000 / eV;
    analysis->histo->ntuData.lastkaonKinEKLDBoost = -1000000 / eV;
    analysis->histo->ntuData.pdgcodeKLDBoost = -1000000;
  }
}

G4bool KLIMAXDegBoostSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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
    if (pname != "") //bkg//
    {
      analysis->histo->ntuData.KLDBooststop[0] = (aStep->GetTrack()->GetPosition())[0] / mm;
      analysis->histo->ntuData.KLDBooststop[1] = (aStep->GetTrack()->GetPosition())[1] / mm;
      analysis->histo->ntuData.KLDBooststop[2] = (aStep->GetTrack()->GetPosition())[2] / mm;

      if (edep > sciEnergyMax) {
        sciEnergyMax = edep;
        TimeKLDBoost = aStep->GetTrack()->GetGlobalTime();
        X = (aStep->GetTrack()->GetPosition())[0] ;
        Y = (aStep->GetTrack()->GetPosition())[1] ;
        Z = (aStep->GetTrack()->GetPosition())[2] ;

        pname = aStep->GetTrack()->GetDefinition()->GetParticleName();
        analysis->histo->ntuData.pdgcodeKLDBoost = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();
        analysis->histo->ntuData.EnergyDepKLDBoost = sciEnergy / eV;
        analysis->histo->ntuData.TimeKLDBoost = TimeKLDBoost / ns;
        analysis->histo->ntuData.XYZKLDBoost[0] = X / mm;
        analysis->histo->ntuData.XYZKLDBoost[1] = Y / mm;
        analysis->histo->ntuData.XYZKLDBoost[2] = Z / mm;
        kaonCounter ++;

        if (kaonCounter == 1)
          analysis->histo->ntuData.kaonKinEKLDBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV; // get the Energy of the first kaon hit

        analysis->histo->ntuData.lastkaonKinEKLDBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV; // get the Energy of the last kaon hit
      }
    }
  }

  return true;
}

void KLIMAXDegBoostSD::EndOfEvent(G4HCofThisEvent*)
{
  kaonCounter = 0;
  kaonKinE = 0;
}
