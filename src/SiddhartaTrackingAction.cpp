#include "../include/SiddhartaDetectorConstruction.h"
#include "../include/SiddhartaTrackingAction.h"
#include "../include/SiddhartaCard.h"

#include <G4TrackingManager.hh>
#include <G4Positron.hh>
#include <G4Track.hh>

SiddhartaTrackingAction::SiddhartaTrackingAction()
{
  G4cout << "-------------------------------------------------------------------------------------------------------" << G4endl;
}

SiddhartaTrackingAction::~SiddhartaTrackingAction()
{
  G4cout << "-------------------------------------------------------------------------------------------------------" << G4endl;
}

void SiddhartaTrackingAction::PreUserTrackingAction(const G4Track* track)
{
  SiddhartaCard* mycard = SiddhartaCard::getInstance();
  int StoreTrajectory = mycard->variables["StoreTrajectory"];

  if (StoreTrajectory != 0)
    fpTrackingManager->SetStoreTrajectory(true);

}

void SiddhartaTrackingAction::PostUserTrackingAction(const G4Track* track) 
{}
