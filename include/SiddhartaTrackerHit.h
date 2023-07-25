#ifndef SiddhartaTrackerHit_h
#define SiddhartaTrackerHit_h 1

#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4Allocator.hh>
#include <G4VHit.hh>

class SiddhartaTrackerHit : public G4VHit
{
public:
  SiddhartaTrackerHit();
  ~SiddhartaTrackerHit();
  SiddhartaTrackerHit(const SiddhartaTrackerHit&);
  const SiddhartaTrackerHit& operator=(const SiddhartaTrackerHit&);
  G4int operator==(const SiddhartaTrackerHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void Draw();
  void Print();

  void SetTrackID(G4int track) {trackID = track;};
  void SetSDDNb(G4int sdd) {SDDNb = sdd;};
  void SetCHIPNb(G4int chip) {CHIPNb = chip;};
  void SetEdep(G4double de) {edep = de;};
  void SetPos(G4ThreeVector xyz) {pos = xyz;};

  G4int GetTrackID() {return trackID;};
  G4int GetSDDNb() {return SDDNb;};
  G4int GetCHIPNb() {return CHIPNb;};
  G4double GetEdep() {return edep;};
  G4ThreeVector GetPos() {return pos;};

private:
  G4int trackID;
  G4int SDDNb;
  G4int CHIPNb;
  G4double edep;
  G4ThreeVector pos;
};

typedef G4THitsCollection<SiddhartaTrackerHit> SiddhartaTrackerHitsCollection;

extern G4Allocator<SiddhartaTrackerHit> SiddhartaTrackerHitAllocator;

inline void* SiddhartaTrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) SiddhartaTrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void SiddhartaTrackerHit::operator delete(void *aHit)
{
  SiddhartaTrackerHitAllocator.FreeSingle((SiddhartaTrackerHit*) aHit);
}

#endif
