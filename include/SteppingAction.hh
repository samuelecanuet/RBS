
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
  

private:
  EventAction *fEventAction;
  RunAction *fRunAction;
  Tracking *fTracking;
  PrimaryGeneratorAction *fGenerator;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
