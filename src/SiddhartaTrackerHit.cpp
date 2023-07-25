#include "../include/SiddhartaTrackerHit.h"

#include <G4VisAttributes.hh>
#include <G4VVisManager.hh>
#include <G4UnitsTable.hh>
#include <G4Colour.hh>
#include <G4Circle.hh>

G4Allocator<SiddhartaTrackerHit> SiddhartaTrackerHitAllocator;

SiddhartaTrackerHit::SiddhartaTrackerHit() {}

SiddhartaTrackerHit::~SiddhartaTrackerHit() {}

SiddhartaTrackerHit::SiddhartaTrackerHit(const SiddhartaTrackerHit& right) : G4VHit()
{
  trackID = right.trackID;
  SDDNb = right.SDDNb;
  edep = right.edep;
  pos = right.pos;
}

const SiddhartaTrackerHit& SiddhartaTrackerHit::operator=(const SiddhartaTrackerHit& right)
{
  trackID = right.trackID;
  SDDNb = right.SDDNb;
  edep = right.edep;
  pos = right.pos;
  return *this;
}

G4int SiddhartaTrackerHit::operator==(const SiddhartaTrackerHit& right) const
{
  return (this==&right) ? 1 : 0;
}

void SiddhartaTrackerHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if (pVVisManager) {
    G4Circle circle(pos);
    circle.SetScreenSize(2.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(0.,1.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void SiddhartaTrackerHit::Print()
{
  G4cout << "  trackID: " << trackID << "  SDDNb: " << SDDNb << "  energy deposit: " << G4BestUnit(edep, "Energy")
         << "  position: " << G4BestUnit(pos, "Length") << G4endl;
}
