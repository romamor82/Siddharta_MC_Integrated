#include "../include/SiddhartaKaonDetectorTopSD.h"
#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaEventAction.h"

#include <G4TrajectoryContainer.hh>
#include <G4EventManager.hh>
#include <G4Trajectory.hh>
#include <G4Event.hh>
#include <G4ios.hh>

SiddhartaEventAction::SiddhartaEventAction() {}

SiddhartaEventAction::~SiddhartaEventAction() {}

void SiddhartaEventAction::BeginOfEventAction(const G4Event*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  analysis->BeginOfEvent();
}

void SiddhartaEventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;

  if (trajectoryContainer)
    n_trajectories = trajectoryContainer->entries();

  if (event_id < 100 || event_id%100 == 0) {
    G4cout << ">>> Event " << evt->GetEventID() << G4endl;
    G4cout << "    " << n_trajectories << " trajectories stored in this event." << G4endl;
  }
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  analysis->EndOfEvent(evt);
}
