#include "../include/SiddhartaLumiDetectorAntiboostSD.h"
#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaHisto.h"

#include <G4HCofThisEvent.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>

SiddhartaLumiDetectorAntiboostSD::SiddhartaLumiDetectorAntiboostSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname = "trackerCollectionLMAntiboost");
}

SiddhartaLumiDetectorAntiboostSD::~SiddhartaLumiDetectorAntiboostSD() {}

void SiddhartaLumiDetectorAntiboostSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new SiddhartaTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
  static G4int HCID = -1;

  if (HCID < 0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection( HCID, trackerCollection );
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    sciEnergy = 0.0;
    sciEnergyMax = 0.0;
    TimeLMAntiboost = -9999.;
    X = -99999.;
    Y = -99999.;
    Z = -99999.;
    analysis->histo->ntuData.particleNameLMAntiboost[0] = '\0';
  }
}

G4bool SiddhartaLumiDetectorAntiboostSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep == 0.)
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
      TimeLMAntiboost = aStep->GetTrack()->GetGlobalTime();
      X = (aStep->GetTrack()->GetPosition())[0] ;
      Y = (aStep->GetTrack()->GetPosition())[1] ;
      Z = (aStep->GetTrack()->GetPosition())[2] ;
      G4String pname ="";
      pname = aStep->GetTrack()->GetDefinition()->GetParticleName();
      int ilen = strlen(pname);

      for (G4int i=0; i<ilen; i++) {
        analysis->histo->ntuData.particleNameLMAntiboost[i] = (pname.data())[i];
      }
      analysis->histo->ntuData.particleNameLMAntiboost[ilen] = '\0';
      analysis->histo->ntuData.TimeLMAntiboost= TimeLMAntiboost/ns;
      analysis->histo->ntuData.XYZLMAntiboost[0] = X/mm;
      analysis->histo->ntuData.XYZLMAntiboost[1] = Y/mm;
      analysis->histo->ntuData.XYZLMAntiboost[2] = Z/mm;
      analysis->histo->ntuData.EnergyDepLMAntiboost= sciEnergy/eV;
    }
  }
  return true;
}

void SiddhartaLumiDetectorAntiboostSD::EndOfEvent(G4HCofThisEvent*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    if (sciEnergy > 0.0)
      analysis->histo->fillHisto("8", sciEnergy, 1.);
    analysis->histo->ntuData.TimeLMAntiboost = TimeLMAntiboost/ns;
    analysis->histo->ntuData.XYZLMAntiboost[0] = X/mm;
    analysis->histo->ntuData.XYZLMAntiboost[1] = Y/mm;
    analysis->histo->ntuData.XYZLMAntiboost[2] = Z/mm;
    analysis->histo->ntuData.EnergyDepLMAntiboost = sciEnergy/eV;
  }
}
