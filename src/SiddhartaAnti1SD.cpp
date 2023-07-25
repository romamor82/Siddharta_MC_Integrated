#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaAnti1SD.h"
#include "../include/SiddhartaHisto.h"

#include <G4HCofThisEvent.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>

SiddhartaAnti1SD::SiddhartaAnti1SD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollectionAnti1");
}

SiddhartaAnti1SD::~SiddhartaAnti1SD() {}

void SiddhartaAnti1SD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new SiddhartaTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
  static G4int HCID = -1;

  if(HCID<0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection(HCID, trackerCollection);
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    for (G4int i=0; i<2; i++) {
      TimeAnti[i] = -9999.;
      EnergyAnti[i] = 0.;
      XYZVertexAnti[i][0] = -9999.;
      XYZVertexAnti[i][1] = -9999.;
      XYZVertexAnti[i][2] = -9999.;
      XYZAnti[i][0] = -9999.;
      XYZAnti[i][1] = -9999.;
      XYZAnti[i][2] = -9999.;
      analysis->histo->ntuData.particleNameAnti[i][0] = '\0';
      analysis->histo->ntuData.parentIDAnti[i] = -1;
      for (unsigned j=0; j<20; j++)
        analysis->histo->ntuData.materialVertexAnti[i][j] = '\0';
      analysis->histo->ntuData.matvertAnti[i] = 0;
    }
  }
}

G4bool SiddhartaAnti1SD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep == 0.)
    return false;

  SiddhartaTrackerHit* newHit = new SiddhartaTrackerHit();
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  trackerCollection->insert(newHit);

  G4int ncopy = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  G4String matname = "";
  matname = aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetMaterial()->GetName();
  int ilen = strlen(matname);

  if (analysis->YesHistos) {
    EnergyAnti[ncopy] += edep;
    if(TimeAnti[ncopy] == -9999.) {
      analysis->histo->ntuData.NbAnti++;

      TimeAnti[ncopy] = aStep->GetTrack()->GetGlobalTime();
      XYZVertexAnti[ncopy][0] = (aStep->GetTrack()->GetVertexPosition())[0];
      XYZVertexAnti[ncopy][1] = (aStep->GetTrack()->GetVertexPosition())[1];
      XYZVertexAnti[ncopy][2] = (aStep->GetTrack()->GetVertexPosition())[2];
      XYZAnti[ncopy][0] = (aStep->GetTrack()->GetPosition())[0];
      XYZAnti[ncopy][1] = (aStep->GetTrack()->GetPosition())[1];
      XYZAnti[ncopy][2] = (aStep->GetTrack()->GetPosition())[2];
      analysis->histo->ntuData.particleAntiPDGEncoding[ncopy] = aStep->GetTrack()->GetDefinition()->GetAntiPDGEncoding();
      analysis->histo->ntuData.parentIDAnti[ncopy] = aStep->GetTrack()->GetParentID();
      analysis->histo->ntuData.KinVertexAnti[ncopy] = (aStep->GetTrack()->GetVertexKineticEnergy()) / eV;

      for (G4int i=0; i<ilen; i++) {
        analysis->histo->ntuData.materialVertexAnti[ncopy][i] = (matname.data())[i];

        if ((i < 11) && ((matname.data())[i] < 123) && ((matname.data())[i] > 31)) {
          analysis->histo->ntuData.matvertAnti[ncopy] = analysis->histo->ntuData.matvertAnti[ncopy]*100 + ((matname.data())[i] - 30);
        }
      }
      analysis->histo->ntuData.parentIDAnti[ncopy] = aStep->GetTrack()->GetParentID();
    }
  }
  return true;
}

void SiddhartaAnti1SD::EndOfEvent(G4HCofThisEvent*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    for(G4int i=0;i<2;i++) {
      analysis->histo->ntuData.TimeAnti[i] = TimeAnti[i] / ns;
      analysis->histo->ntuData.EnergyAnti[i] = EnergyAnti[i] / eV;
      analysis->histo->ntuData.XYZVertexAnti[i][0] = XYZVertexAnti[i][0] / mm;
      analysis->histo->ntuData.XYZVertexAnti[i][1] = XYZVertexAnti[i][1] / mm;
      analysis->histo->ntuData.XYZVertexAnti[i][2] = XYZVertexAnti[i][2] / mm;
      analysis->histo->ntuData.XYZAnti[i][0] = XYZAnti[i][0] / mm;
      analysis->histo->ntuData.XYZAnti[i][1] = XYZAnti[i][1] / mm;
      analysis->histo->ntuData.XYZAnti[i][2] = XYZAnti[i][2] / mm;
    }
  }
}
