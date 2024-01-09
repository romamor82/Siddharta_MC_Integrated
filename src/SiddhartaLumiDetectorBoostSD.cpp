#include "../include/SiddhartaLumiDetectorBoostSD.h"
#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaHisto.h"

#include <G4HCofThisEvent.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>

SiddhartaLumiDetectorBoostSD::SiddhartaLumiDetectorBoostSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname = "trackerCollectionLMBoost");
}

SiddhartaLumiDetectorBoostSD::~SiddhartaLumiDetectorBoostSD() {}

void SiddhartaLumiDetectorBoostSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new SiddhartaTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
  static G4int HCID = -1;

  if (HCID < 0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection(HCID, trackerCollection);
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    sciEnergy = 0.0;
    sciEnergyMax = 0.0;
    TimeLMBoost = -9999.;
    X = -99999. ;
    Y = -99999. ;
    Z = -99999. ;
    analysis->histo->ntuData.particleNameLMBoost[0] = '\0';
    analysis->histo->ntuData.pdgcodeLMBoost = -1000000;
    analysis->histo->ntuData.EnergyDepLMBoost = -1000000./eV;
    analysis->histo->ntuData.TimeLMBoost = -1000000./ns;
    analysis->histo->ntuData.XYZLMBoost[0] = -1000000./mm;
    analysis->histo->ntuData.XYZLMBoost[1] = -1000000./mm;
    analysis->histo->ntuData.XYZLMBoost[2] = -1000000./mm;
    analysis->histo->ntuData.XYZLMBoostKaonstop[0] = -1000000./mm;
    analysis->histo->ntuData.XYZLMBoostKaonstop[1] = -1000000./mm;
    analysis->histo->ntuData.XYZLMBoostKaonstop[2] = -1000000./mm;
    analysis->histo->ntuData.lastkaonKinELMBoost = -1000000./eV;
  }
}

G4bool SiddhartaLumiDetectorBoostSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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
      TimeLMBoost = aStep->GetTrack()->GetGlobalTime();
      X = (aStep->GetTrack()->GetPosition())[0] ;
      Y = (aStep->GetTrack()->GetPosition())[1] ;
      Z = (aStep->GetTrack()->GetPosition())[2] ;
      G4String pname ="";
      pname = aStep->GetTrack()->GetDefinition()->GetParticleName();
      int ilen = strlen(pname);

      for (G4int i=0; i<ilen; i++) {
        analysis->histo->ntuData.particleNameLMBoost[i] = (pname.data())[i];
      }
      analysis->histo->ntuData.pdgcodeLMBoost = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();
      analysis->histo->ntuData.particleNameLMBoost[ilen] = '\0';
      analysis->histo->ntuData.EnergyDepLMBoost = sciEnergy/eV;
      analysis->histo->ntuData.TimeLMBoost = TimeLMBoost/ns;
      analysis->histo->ntuData.XYZLMBoost[0] = X/mm;
      analysis->histo->ntuData.XYZLMBoost[1] = Y/mm;
      analysis->histo->ntuData.XYZLMBoost[2] = Z/mm;
    }
    if (aStep->GetTrack()->GetDynamicParticle()->GetPDGcode() == -321) // KAONS
    {
      kaonCounter++;

      if (kaonCounter == 1) {analysis->histo->ntuData.kaonKinELMBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV;} // get the Energy of the first kaon hit

      analysis->histo->ntuData.lastkaonKinELMBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV; // get the Energy of the last kaon hit
      analysis->histo->ntuData.XYZLMBoostKaonstop[0] = X/mm;
      analysis->histo->ntuData.XYZLMBoostKaonstop[1] = Y/mm;
      analysis->histo->ntuData.XYZLMBoostKaonstop[2] = Z/mm;
    }
  }

  return true;

  }

void SiddhartaLumiDetectorBoostSD::EndOfEvent(G4HCofThisEvent*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  if (analysis->YesHistos) {
    if (sciEnergy > 0.0)
      analysis->histo->fillHisto("7", sciEnergy, 1.);
    analysis->histo->ntuData.EnergyDepLMBoost = sciEnergy/eV;
    analysis->histo->ntuData.TimeLMBoost = TimeLMBoost/ns;
    analysis->histo->ntuData.XYZLMBoost[0] = X/mm;
    analysis->histo->ntuData.XYZLMBoost[1] = Y/mm;
    analysis->histo->ntuData.XYZLMBoost[2] = Z/mm;
  }
  kaonCounter = 0;
}
