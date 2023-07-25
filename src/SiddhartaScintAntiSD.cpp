#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaScintAntiSD.h"
#include "../include/SiddhartaHisto.h"

#include <G4HCofThisEvent.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>

SiddhartaScintAntiSD::SiddhartaScintAntiSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname = "trackerCollectionScintAnti");
}

SiddhartaScintAntiSD::~SiddhartaScintAntiSD() {}

void SiddhartaScintAntiSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new SiddhartaTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
  static G4int HCID = -1;

  if (HCID < 0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection( HCID, trackerCollection );
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    nHitFull = 0;
    for (G4int i=0; i<48; i++) {
      nHit[i] = 0;
      EnergyDep[i] = 0.;
      TimeScintAnti[i] = -9999.;
      Copy[i] = -1;	//Copy following numeration by nHitFull!!!
      analysis->histo->ntuData.particleScintAntiPDGEncoding[i] = 0;
      analysis->histo->ntuData.matVertIDFromMap_Veto2[i] = 0;
    }
  }
}

G4bool SiddhartaScintAntiSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep == 0.)
    return false;

  G4String matname = "";
  matname = aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetMaterial()->GetName();

  SiddhartaTrackerHit* newHit = new SiddhartaTrackerHit();
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  trackerCollection->insert(newHit);
  G4int ncopy = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1)->GetCopyNo();
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  if (analysis->YesHistos) {
    EnergyDep[nHitFull] += edep;	//Time and energy following numeration by Copy!!!
//std::cout << ncopy << " " << edep << " " << nHitFull << std::endl;
    nHit[ncopy-1] += 1;
    if (nHit[ncopy-1] == 1) {
      if(TimeScintAnti[nHitFull] == -9999.) {
        TimeScintAnti[nHitFull] = aStep->GetTrack()->GetGlobalTime();	//Time and energy following numeration by Copy!!!
        Copy[nHitFull] = ncopy;	//Copy following numeration by nHitFull!!!
        analysis->histo->ntuData.particleScintAntiPDGEncoding[nHitFull] = aStep->GetTrack()->GetDefinition()->GetAntiPDGEncoding();
        analysis->histo->ntuData.matVertIDFromMap_Veto2[nHitFull] = analysis->histo->getMatIDFromString(matname);
//std::cout << aStep->GetTrack()->GetGlobalTime()/ns << " " << matname << std::endl;
        nHitFull += 1;
      }
    }
  }
  return true;
}

void SiddhartaScintAntiSD::EndOfEvent(G4HCofThisEvent*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    analysis->histo->ntuData.nHitSciAnti = nHitFull;
    for (G4int i=0; i<nHitFull; i++) {
      analysis->histo->ntuData.TimeScintAnti[i] = TimeScintAnti[i]/ns;	//Time and energy following numeration by Copy!!!
      analysis->histo->ntuData.EnergyScintAnti[i] = EnergyDep[i]/eV;	//Time and energy following numeration by Copy!!!
      analysis->histo->ntuData.copyScintAnti[i] = Copy[i];	//Copy following numeration by nHitFull!!!
      analysis->histo->fillHistogram("Veto2_time", TimeScintAnti[i]/ns/ns, doubleCheck(0), doubleCheck(0));
/*std::cout << TimeScintAnti[i]/ns << " " << EnergyDep[i]/eV << " " << Copy[i] << std::endl;
int n;
std::cin >> n;*/
    }
  }
}
