#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaRunAction.h"
#include "../include/SiddhartaCard.h"

#include <G4Run.hh>

#include <ctime>
#include <map>

SiddhartaRunAction::SiddhartaRunAction() {}

SiddhartaRunAction::~SiddhartaRunAction() {}

void SiddhartaRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  SiddhartaCard* mycard = SiddhartaCard::getInstance();

  time_t t=time(0);
  int seed = t - 1299000000;
  mycard->SetG4RandGaussSeed(seed);

  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  analysis->BeginOfRun();

  G4int RunN = aRun->GetRunID();
  if (RunN % 1000 == 0)
  G4cout << "### Run : " << RunN << G4endl;
}

void SiddhartaRunAction::EndOfRunAction(const G4Run*)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  analysis->EndOfRun();
}
