#ifndef TrackingAction_h
#define TrackingAction_h 1

#include <G4UserTrackingAction.hh>
#include <globals.hh>

class SiddhartaTrackingAction : public G4UserTrackingAction
{
public:
  SiddhartaTrackingAction();
  ~SiddhartaTrackingAction();

  void  PreUserTrackingAction(const G4Track*);
  void PostUserTrackingAction(const G4Track*);
};

#endif
