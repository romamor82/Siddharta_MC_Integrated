#include "../include/SiddhartaAnalysisManager.h"
#include "../include/KLIMAXDegAntiBoostSD.h"
#include "../include/SiddhartaHisto.h"

#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>


KLIMAXDegAntiBoostSD::KLIMAXDegAntiBoostSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollectionKLDAntiBoost");
}

KLIMAXDegAntiBoostSD::~KLIMAXDegAntiBoostSD() {}

void KLIMAXDegAntiBoostSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new SiddhartaTrackerHitsCollection(SensitiveDetectorName,collectionName[0]);
  static G4int HCID = -1;

  if (HCID<0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection( HCID, trackerCollection ); 

  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    sciEnergy = 0.0;
    sciEnergyMax = 0.0;
    TimeKLDAntiBoost = -9999.;
    X = -99999.;
    Y = -99999.;
    Z = -99999.;
    analysis->histo->ntuData.KLDAntiBooststop[0] = -1000000 / mm;
    analysis->histo->ntuData.KLDAntiBooststop[1] = -1000000 /  mm;
    analysis->histo->ntuData.KLDAntiBooststop[2] = -1000000 / mm;
    analysis->histo->ntuData.EnergyDepKLDAntiBoost = -1000000 / eV;
    analysis->histo->ntuData.TimeKLDAntiBoost = -1000000 / ns;
    analysis->histo->ntuData.XYZKLDAntiBoost[0] = -1000000 / mm;
    analysis->histo->ntuData.XYZKLDAntiBoost[1] = -1000000 / mm;
    analysis->histo->ntuData.XYZKLDAntiBoost[2] = -1000000 / mm;
    analysis->histo->ntuData.kaonKinEKLDAntiBoost = -1000000 / eV;
    analysis->histo->ntuData.lastkaonKinEKLDAntiBoost = -1000000 / eV;
    analysis->histo->ntuData.pdgcodeKLDAntiBoost = -1000000;
  }
}

G4bool KLIMAXDegAntiBoostSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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
      analysis->histo->ntuData.KLDAntiBooststop[0] = (aStep->GetTrack()->GetPosition())[0] / mm;
      analysis->histo->ntuData.KLDAntiBooststop[1] = (aStep->GetTrack()->GetPosition())[1] / mm;
      analysis->histo->ntuData.KLDAntiBooststop[2] = (aStep->GetTrack()->GetPosition())[2] / mm;

      if (edep > sciEnergyMax) {
        sciEnergyMax = edep;
        TimeKLDAntiBoost = aStep->GetTrack()->GetGlobalTime();
        X = (aStep->GetTrack()->GetPosition())[0];
        Y = (aStep->GetTrack()->GetPosition())[1];
        Z = (aStep->GetTrack()->GetPosition())[2];

        pname = aStep->GetTrack()->GetDefinition()->GetParticleName();
        analysis->histo->ntuData.pdgcodeKLDAntiBoost = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();
        analysis->histo->ntuData.EnergyDepKLDAntiBoost = sciEnergy / eV;
        analysis->histo->ntuData.TimeKLDAntiBoost = TimeKLDAntiBoost / ns;
        analysis->histo->ntuData.XYZKLDAntiBoost[0] = X / mm;
        analysis->histo->ntuData.XYZKLDAntiBoost[1] = Y / mm;
        analysis->histo->ntuData.XYZKLDAntiBoost[2] = Z / mm;
        kaonCounter ++;

        if (kaonCounter == 1)
          analysis->histo->ntuData.kaonKinEKLDAntiBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV; // get the Energy of the first kaon hit

        analysis->histo->ntuData.lastkaonKinEKLDAntiBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV; // get the Energy of the last kaon hit
      }
    }
  }

  return true;
}


void KLIMAXDegAntiBoostSD::EndOfEvent(G4HCofThisEvent*)
{
  kaonCounter = 0;
  kaonKinE = 0;
}
