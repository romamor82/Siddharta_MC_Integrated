//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

#ifndef SiddhartaPhysicsList_h
#define SiddhartaPhysicsList_h 1

#include <G4VUserPhysicsList.hh>
#include <G4SystemOfUnits.hh>
#include <globals.hh>

class SiddhartaPhysicsList: public G4VUserPhysicsList
{
public:
  SiddhartaPhysicsList();
  ~SiddhartaPhysicsList();

  void SetCuts();
  void SetGammaCut(G4double);
  void SetElectronCut(G4double);
  void SetPositronCut(G4double);

  void SetGammaLowLimit(G4double);
  void SetElectronLowLimit(G4double);

  void SetGELowLimit(G4double);

  void SetLowEnSecPhotCut(G4double);
  void SetLowEnSecElecCut(G4double);

  void SetProtonCut(G4double);
  void SetNeutronCut(G4double);

private:
  G4double LowEnCut;
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  G4double cutForKaonPlus;
  G4double cutForKaonMinus;
  G4double cutForPionPlus;
  G4double cutForPionMinus;
  G4double cutForProton;
  G4double cutForNeutron;

protected:
  void ConstructParticle();
  void ConstructProcess();

  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();
  void ConstructShortLived();
  void ConstructIons();

  void AddTransportation();
  void ConstructGeneral();
  void ConstructEM();
  void ConstructHadronic();
  void AddStepMax();
};

#endif
