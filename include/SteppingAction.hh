
/// \file SteppingAction.hh
/// \brief Definition of the SteppingAction class

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "TH1D.h"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4EmCalculator.hh"
#include "G4Exception.hh"


class EventAction;
class RunAction;
class Tracking;
class PrimaryGeneratorAction;
class Messenger;

/// Stepping action class
///

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(EventAction *eventAction, RunAction *runAction, Tracking *tracking, PrimaryGeneratorAction *generator);
  virtual ~SteppingAction();

  // method from the base class
  virtual void UserSteppingAction(const G4Step *);
  G4double CollisionEnergy(G4double Energy, const G4Element *Element, G4double Angle);
  G4ThreeVector DirectionRBS();
  G4double CrossSectionRBS(G4double Energy, const G4Element *Element);
  const G4Element *RandomElement(G4Material *Material);

  TH1D *Hist_sin4 = new TH1D("sin4", "sin4", 10 * pi, (135 - 7 - 2 * 5) * CLHEP::deg, (135 + 7 + 2 * 5) * CLHEP::deg);


private:
  EventAction *fEventAction;
  RunAction *fRunAction;
  Tracking *fTracking;
  PrimaryGeneratorAction *fGenerator;

  Messenger* msg;
  
  G4EmCalculator emCalculator;
  G4double energy = 0;
  G4int Parent;
  G4ParticleDefinition* InitialParticle;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
