#ifndef SiddhartaRunAction_h
#define SiddhartaRunAction_h 1

#include "SiddhartaAnalysisManager.h"
#include "SiddhartaHisto.h"

#include <G4UserRunAction.hh>
#include <globals.hh>

class G4Run;

class SiddhartaRunAction : public G4UserRunAction
{
public:
  SiddhartaRunAction();
  ~SiddhartaRunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

private:
  SiddhartaHisto* histosManager;
};

#endif
