#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaAnti3SD.h"
#include "../include/SiddhartaHisto.h"
#include "../include/SiddhartaCard.h" //??

#include <G4HCofThisEvent.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>

SiddhartaAnti3SD::SiddhartaAnti3SD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollectionAnti3");
}

SiddhartaAnti3SD::~SiddhartaAnti3SD() {}

void SiddhartaAnti3SD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new SiddhartaTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
  static G4int HCID = -1;

  if (HCID<0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection( HCID, trackerCollection );

  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  SiddhartaCard* mycard = SiddhartaCard::getInstance();
  int SiddhartaSetup = mycard->variables["SiddhartaSetupVersion"];
  int ntot = 2;

  if (SiddhartaSetup == 6)
    ntot = 2;
  else if (SiddhartaSetup == 7 || SiddhartaSetup == 8)
    ntot = 12;

  if (analysis->YesHistos) {
    for (G4int i=0; i<12; i++) {
      TimeAnti[i] = -9999.;
      EnergyAnti[i] = 0.;
      XYZVertexAnti[i][0] = -9999.;
      XYZVertexAnti[i][1] = -9999.;
      XYZVertexAnti[i][2] = -9999.;
      XYZAnti[i][0] = -9999.;
      XYZAnti[i][1] = -9999.;
      XYZAnti[i][2] = -9999.;
      analysis->histo->ntuData.NbAnti = 0;
      analysis->histo->ntuData.particleNameAnti[5+i][0] = '\0';//[2+i]
      analysis->histo->ntuData.parentIDAnti[5+i] = -1;//[2+i]
      for (unsigned j=0; j<20; j++)
        analysis->histo->ntuData.materialVertexAnti[5+i][j] = '\0';
      analysis->histo->ntuData.matvertAnti[5+i] = 0;
      analysis->histo->ntuData.matVertIDFromMap_Veto1[5+i] = 0;
    }
  }

}

G4bool SiddhartaAnti3SD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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
      analysis->histo->ntuData.particleAntiPDGEncoding[5+ncopy] = aStep->GetTrack()->GetDefinition()->GetAntiPDGEncoding();
      analysis->histo->ntuData.parentIDAnti[5+ncopy] = aStep->GetTrack()->GetParentID();
      analysis->histo->ntuData.KinVertexAnti[5+ncopy] = (aStep->GetTrack()->GetVertexKineticEnergy()) / eV;

      for (G4int i=0; i<ilen; i++) {
        analysis->histo->ntuData.materialVertexAnti[5+ncopy][i] = (matname.data())[i];

        if ((i < 11) && ((matname.data())[i] < 123) && ((matname.data())[i] > 31)) {
          analysis->histo->ntuData.matvertAnti[5+ncopy] = analysis->histo->ntuData.matvertAnti[5+ncopy]*100 + ((matname.data())[i] - 30);
        }
      }
      analysis->histo->ntuData.matVertIDFromMap_Veto1[5+ncopy] = analysis->histo->getMatIDFromString(matname);
    }
  }
  return true;
}

void SiddhartaAnti3SD::EndOfEvent(G4HCofThisEvent*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  SiddhartaCard* mycard = SiddhartaCard::getInstance();
  int SiddhartaSetup = mycard->variables["SiddhartaSetupVersion"];

  if (SiddhartaSetup == 7 || SiddhartaSetup == 8 || SiddhartaSetup == 2020 || SiddhartaSetup == 2023) {
    if (analysis->YesHistos) {
      for(G4int i=0;i<12;i++) {
        analysis->histo->ntuData.TimeAnti[5+i] = TimeAnti[i] / ns;
        analysis->histo->ntuData.EnergyAnti[5+i] = EnergyAnti[i] / eV;
        analysis->histo->ntuData.XYZVertexAnti[5+i][0] = XYZVertexAnti[i][0] / mm;
        analysis->histo->ntuData.XYZVertexAnti[5+i][1] = XYZVertexAnti[i][1] / mm;
        analysis->histo->ntuData.XYZVertexAnti[5+i][2] = XYZVertexAnti[i][2] / mm;
        analysis->histo->ntuData.XYZAnti[5+i][0] = XYZAnti[i][0] / mm;
        analysis->histo->ntuData.XYZAnti[5+i][1] = XYZAnti[i][1] / mm;
        analysis->histo->ntuData.XYZAnti[5+i][2] = XYZAnti[i][2] / mm;
      }
    }
  } else if (analysis->YesHistos) {
    for(G4int i=0;i<2;i++) {
      analysis->histo->ntuData.TimeAnti[5+i] = TimeAnti[i] / ns;
      analysis->histo->ntuData.EnergyAnti[5+i] = EnergyAnti[i] / eV;
      analysis->histo->ntuData.XYZVertexAnti[5+i][0] = XYZVertexAnti[i][0] / mm;
      analysis->histo->ntuData.XYZVertexAnti[5+i][1] = XYZVertexAnti[i][1] / mm;
      analysis->histo->ntuData.XYZVertexAnti[5+i][2] = XYZVertexAnti[i][2] / mm;
      analysis->histo->ntuData.XYZAnti[5+i][0] = XYZAnti[i][0] / mm;
      analysis->histo->ntuData.XYZAnti[5+i][1] = XYZAnti[i][1] / mm;
      analysis->histo->ntuData.XYZAnti[5+i][2] = XYZAnti[i][2] / mm;
    }
  }
}
