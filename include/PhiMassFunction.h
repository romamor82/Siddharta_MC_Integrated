///////////////////////////
//  PhiMassFunction.h
//   2021.11.23
//    H. Shi
///////////////////////////

#ifndef _PHIMASSFUNCTION_HH
#define _PHIMASSFUNCTION_HH

#include <G4SystemOfUnits.hh>
#include <globals.hh>

#include <TMath.h>

const static G4double mphi  = 1019.46*MeV;
const static G4double phiGm = 4.26*MeV;

G4double phiMassFunc(G4double*, G4double*);

#endif
