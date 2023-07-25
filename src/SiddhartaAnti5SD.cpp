#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaAnti5SD.h"
#include "../include/SiddhartaHisto.h"

#include <G4HCofThisEvent.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>

SiddhartaAnti5SD::SiddhartaAnti5SD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollectionAnti5");
}

SiddhartaAnti5SD::~SiddhartaAnti5SD(){ }

void SiddhartaAnti5SD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new SiddhartaTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
  static G4int HCID = -1;

  if (HCID<0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection( HCID, trackerCollection );
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    TimeAnti = -9999.;
    EnergyAnti = 0.;
    XYZVertexAnti[0] = -9999.;
    XYZVertexAnti[1] = -9999.;
    XYZVertexAnti[2] = -9999.;
    XYZAnti[0] = -9999.;
    XYZAnti[1] = -9999.;
    XYZAnti[2] = -9999.;
    analysis->histo->ntuData.particleNameAnti[4][0] = '\0';//[8]
    analysis->histo->ntuData.parentIDAnti[4] = -1;//[8]
    for (unsigned j=0; j<20; j++)
      analysis->histo->ntuData.materialVertexAnti[4][j] = '\0';
    analysis->histo->ntuData.matvertAnti[4] = 0;
  }
}

G4bool SiddhartaAnti5SD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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

  G4String matname = "";
  matname = aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetMaterial()->GetName();
  int ilen = strlen(matname);

  if (analysis->YesHistos) {
    EnergyAnti += edep;
    if (TimeAnti == -9999.) {
      TimeAnti = aStep->GetTrack()->GetGlobalTime();
      XYZVertexAnti[0] = (aStep->GetTrack()->GetVertexPosition())[0];
      XYZVertexAnti[1] = (aStep->GetTrack()->GetVertexPosition())[1];
      XYZVertexAnti[2] = (aStep->GetTrack()->GetVertexPosition())[2];
      XYZAnti[0] = (aStep->GetTrack()->GetPosition())[0];
      XYZAnti[1] = (aStep->GetTrack()->GetPosition())[1];
      XYZAnti[2] = (aStep->GetTrack()->GetPosition())[2];
      analysis->histo->ntuData.particleAntiPDGEncoding[4] = aStep->GetTrack()->GetDefinition()->GetAntiPDGEncoding();
      analysis->histo->ntuData.parentIDAnti[4] = aStep->GetTrack()->GetParentID();
      analysis->histo->ntuData.KinVertexAnti[4] = (aStep->GetTrack()->GetVertexKineticEnergy()) / eV;

      for (G4int i=0; i<ilen; i++) {
        analysis->histo->ntuData.materialVertexAnti[4][i] = (matname.data())[i];

        if ((i < 11) && ((matname.data())[i] < 123) && ((matname.data())[i] > 31)) {
          analysis->histo->ntuData.matvertAnti[4] = analysis->histo->ntuData.matvertAnti[4]*100 + ((matname.data())[i] - 30);
        }
      }
    }
  }
  return true;
}

void SiddhartaAnti5SD::EndOfEvent(G4HCofThisEvent*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    analysis->histo->ntuData.TimeAnti[4] = TimeAnti / ns;
    analysis->histo->ntuData.EnergyAnti[4] = EnergyAnti / eV;
    analysis->histo->ntuData.XYZVertexAnti[4][0] = XYZVertexAnti[0] / mm;
    analysis->histo->ntuData.XYZVertexAnti[4][1] = XYZVertexAnti[1] / mm;
    analysis->histo->ntuData.XYZVertexAnti[4][2] = XYZVertexAnti[2] / mm;
    analysis->histo->ntuData.XYZAnti[4][0] = XYZAnti[0] / mm;
    analysis->histo->ntuData.XYZAnti[4][1] = XYZAnti[1] / mm;
    analysis->histo->ntuData.XYZAnti[4][2] = XYZAnti[2] / mm;
  }
}
