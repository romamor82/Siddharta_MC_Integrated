#include "../include/SiddhartaLumiDetectorAntiBoostSD.h"
#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaHisto.h"

#include <G4HCofThisEvent.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4VSolid.hh>
#include <G4Step.hh>
#include <G4ios.hh>

SiddhartaLumiDetectorAntiBoostSD::SiddhartaLumiDetectorAntiBoostSD(G4String name) : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname = "trackerCollectionLMAntiBoost");
}

SiddhartaLumiDetectorAntiBoostSD::~SiddhartaLumiDetectorAntiBoostSD() {}

void SiddhartaLumiDetectorAntiBoostSD::Initialize(G4HCofThisEvent* HCE)
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
    TimeLMAntiBoost = -9999.;
    X = -99999.;
    Y = -99999.;
    Z = -99999.;
    analysis->histo->ntuData.particleNameLMAntiBoost[0] = '\0';
    analysis->histo->ntuData.pdgcodeLMAntiBoost = -1000000;
    analysis->histo->ntuData.EnergyDepLMAntiBoost = -1000000./eV;
    analysis->histo->ntuData.TimeLMAntiBoost = -1000000./ns;
    analysis->histo->ntuData.XYZLMAntiBoost[0] = -1000000./mm;
    analysis->histo->ntuData.XYZLMAntiBoost[1] = -1000000./mm;
    analysis->histo->ntuData.XYZLMAntiBoost[2] = -1000000./mm;
    analysis->histo->ntuData.XYZLMAntiBoostKaonstop[0] = -1000000./mm;
    analysis->histo->ntuData.XYZLMAntiBoostKaonstop[1] = -1000000./mm;
    analysis->histo->ntuData.XYZLMAntiBoostKaonstop[2] = -1000000./mm;
    analysis->histo->ntuData.lastkaonKinELMAntiBoost = -1000000./eV;
  }
}

G4bool SiddhartaLumiDetectorAntiBoostSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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
	G4String pname ="";
	pname = aStep->GetTrack()->GetDefinition()->GetParticleName();

	sciEnergy += edep;

	if (analysis->YesHistos) {
		X = (aStep->GetTrack()->GetPosition())[0] ;
		Y = (aStep->GetTrack()->GetPosition())[1] ;
		Z = (aStep->GetTrack()->GetPosition())[2] ;
		if (edep > sciEnergyMax) {
			sciEnergyMax = edep;
			TimeLMAntiBoost = aStep->GetTrack()->GetGlobalTime();

			analysis->histo->ntuData.pdgcodeLMAntiBoost = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();
			analysis->histo->ntuData.XYZLMAntiBoost[0] = X/mm;
			analysis->histo->ntuData.XYZLMAntiBoost[1] = Y/mm;
			analysis->histo->ntuData.XYZLMAntiBoost[2] = Z/mm;
			analysis->histo->ntuData.EnergyDepLMAntiBoost= sciEnergy/eV;
			analysis->histo->ntuData.TimeLMAntiBoost = TimeLMAntiBoost/ns;
		}
		if(aStep->GetTrack()->GetDynamicParticle()->GetPDGcode()==-321) // KAONS
		{
//			G4cout << "K- in LumiAB " << (aStep->GetTrack()->GetKineticEnergy()) / eV << G4endl;
			kaonCounter++;
			if(kaonCounter == 1) {analysis->histo->ntuData.kaonKinELMAntiBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV;} // get the Energy of the first kaon hit
			analysis->histo->ntuData.lastkaonKinELMAntiBoost = (aStep->GetTrack()->GetKineticEnergy()) / eV; // get the Energy of the last kaon hit
			analysis->histo->ntuData.XYZLMAntiBoostKaonstop[0] = X/mm;
			analysis->histo->ntuData.XYZLMAntiBoostKaonstop[1] = Y/mm;
			analysis->histo->ntuData.XYZLMAntiBoostKaonstop[2] = Z/mm;
		}
	}
	return true;
}

void SiddhartaLumiDetectorAntiBoostSD::EndOfEvent(G4HCofThisEvent*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();

  if (analysis->YesHistos) {
    if (sciEnergy > 0.0)
      analysis->histo->fillHisto("8", sciEnergy, 1.);
    analysis->histo->ntuData.TimeLMAntiBoost = TimeLMAntiBoost/ns;
    analysis->histo->ntuData.XYZLMAntiBoost[0] = X/mm;
    analysis->histo->ntuData.XYZLMAntiBoost[1] = Y/mm;
    analysis->histo->ntuData.XYZLMAntiBoost[2] = Z/mm;
    analysis->histo->ntuData.EnergyDepLMAntiBoost = sciEnergy/eV;
  }
  kaonCounter = 0;
}
