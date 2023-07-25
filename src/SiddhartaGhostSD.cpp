#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaGhostSD.h"
#include "../include/SiddhartaHisto.h"

#include <G4HCofThisEvent.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>

SiddhartaGhostSD::SiddhartaGhostSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollectionGhost");
}

SiddhartaGhostSD::~SiddhartaGhostSD() {}

void SiddhartaGhostSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new SiddhartaTrackerHitsCollection(SensitiveDetectorName,collectionName[0]);
  static G4int HCID = -1;
  if(HCID<0)
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

  HCE->AddHitsCollection( HCID, trackerCollection );
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    TimeGhost = -9999.;
    X = -99999.;
    Y = -99999.;
    Z = -99999.;
    sciEnergy = 0.0;
    sciEnergyMax = 0.0;
    KineticGhost = 0.0;
  }
}

G4bool SiddhartaGhostSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep==0.)
    return false;

  SiddhartaTrackerHit* newHit = new SiddhartaTrackerHit();
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  trackerCollection->insert(newHit);
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  sciEnergy += edep;

  if (analysis->YesHistos) {
    if (edep > sciEnergyMax) {
      sciEnergyMax = edep;
      TimeGhost = aStep->GetTrack()->GetGlobalTime();
      KineticGhost = aStep->GetTrack()->GetKineticEnergy();
      X = (aStep->GetTrack()->GetPosition())[0] ;
      Y = (aStep->GetTrack()->GetPosition())[1] ;
      Z = (aStep->GetTrack()->GetPosition())[2] ;
    }
  }
  return true;
}

void SiddhartaGhostSD::EndOfEvent(G4HCofThisEvent*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    analysis->histo->ntuData.TimeGhost = TimeGhost/ns;
    analysis->histo->ntuData.KineticGhost = KineticGhost/eV;
    analysis->histo->ntuData.XYZGhost[0] = X/mm;
    analysis->histo->ntuData.XYZGhost[1] = Y/mm;
    analysis->histo->ntuData.XYZGhost[2] = Z/mm;
    analysis->histo->ntuData.EnergyDepGhost = sciEnergy/eV;
  }
}
