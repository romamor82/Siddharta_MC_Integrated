#ifndef SiddhartaEventAction_h
#define SiddhartaEventAction_h 1

#include <G4UserEventAction.hh>

class G4Event;

class SiddhartaEventAction : public G4UserEventAction
{
public:
  SiddhartaEventAction();
 ~SiddhartaEventAction();

public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
};

#endif
