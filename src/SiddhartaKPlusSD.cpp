#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaKPlusSD.h"
#include "../include/SiddhartaHisto.h"

#include <G4HCofThisEvent.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>

SiddhartaKPlusSD::SiddhartaKPlusSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname = "trackerCollectionKPlus");
}

SiddhartaKPlusSD::~SiddhartaKPlusSD() {}

void SiddhartaKPlusSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new SiddhartaTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
  static G4int HCID = -1;

  if (HCID < 0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection(HCID, trackerCollection );
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    sciEnergy = 0.0;
    sciEnergyMax = 0.0;
    TimeKPlusDetector = -9999.;
    X = -99999.;
    Y = -99999.;
    Z = -99999.;
    analysis->histo->ntuData.particleNameKPlusDetector[0] = '\0';
  }
}

G4bool SiddhartaKPlusSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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
      TimeKPlusDetector = aStep->GetTrack()->GetGlobalTime();
      X = (aStep->GetTrack()->GetPosition())[0] ;
      Y = (aStep->GetTrack()->GetPosition())[1] ;
      Z = (aStep->GetTrack()->GetPosition())[2] ;
      G4String pname ="";
      pname = aStep->GetTrack()->GetDefinition()->GetParticleName();
      int ilen = strlen(pname);

      for (G4int i=0; i<ilen; i++) {
        analysis->histo->ntuData.particleNameKPlusDetector[i] = (pname.data())[i];
      }
      analysis->histo->ntuData.particleNameKPlusDetector[ilen] = '\0';
      analysis->histo->ntuData.TimeKPlusDetector = TimeKPlusDetector/ns;
      analysis->histo->ntuData.XYZKPlusDetector[0] = X/mm;
      analysis->histo->ntuData.XYZKPlusDetector[1] = Y/mm;
      analysis->histo->ntuData.XYZKPlusDetector[2] = Z/mm;
      analysis->histo->ntuData.EnergyDepKPlusDetector = sciEnergy/eV;
    }
  }
  return true;
}

void SiddhartaKPlusSD::EndOfEvent(G4HCofThisEvent*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  if (analysis->YesHistos) {
    if (sciEnergy > 0.0)
      analysis->histo->fillHistogram("KaonPlusDetector_time", TimeKPlusDetector/ns, doubleCheck(0), doubleCheck(0));

    analysis->histo->ntuData.TimeKPlusDetector = TimeKPlusDetector/ns;
    analysis->histo->ntuData.XYZKPlusDetector[0] = X/mm;
    analysis->histo->ntuData.XYZKPlusDetector[1] = Y/mm;
    analysis->histo->ntuData.XYZKPlusDetector[2] = Z/mm;
    analysis->histo->ntuData.EnergyDepKPlusDetector = sciEnergy/eV;
  }
}
