///////////////////////////////
//   PhiMassFunction.cc
//     2021.11.23
//     H. Shi
//////////////////////////////

#include  "../include/PhiMassFunction.h"

/*
Unneccessary file -> Move to more general one
 */

// Phi Mass PDF
G4double phiMassFunc( G4double *x, G4double *par )
{
  G4double amp;
  amp = par[1]/TMath::Pi()/((x[0]-par[0])*(x[0]-par[0])+par[1]*par[1]);
  par[0] = mphi;
  par[1] = phiGm/2;

  return amp;
}
