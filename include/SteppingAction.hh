
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
  G4double CrossSectionRBS(G4double Energy, const G4Element *Element, G4double theta);
  const G4Element *RandomElement(G4Material *Material);

  TH1D *Hist_sin4 = new TH1D("sin4", "sin4", 10 * pi, (135 - 7 - 2 * 5) * CLHEP::deg, (135 + 7 + 2 * 5) * CLHEP::deg);



//////////////::::::::::
G4double CalcAngleCMFrame(G4double angle, G4double M1, G4double M2);
	//function to calculate energy in the CM reference frame
    G4double CalcEnergyCMFrame(G4double energy, G4double M1, G4double M2);	       
	// RBS xsec in the Lab frame from the CM reference fram
    G4double CalcDiffRuthXsecLAB(G4double M1, G4double M2, G4double angle, G4double xsection);
        G4double CalcDiffRuthXsecCM(G4double E, G4double angleCM ,G4double Z1, G4double Z2); 

    	// function for Total RBS yield, combining other functions into single one
    G4double CalculateTotalRBSYield(G4double energy, G4double M1, G4double M2, G4double Z1, G4double Z2, G4double angle);
      G4double CalcAndersenScreening(G4double energy_cm, G4double angle_cm, G4double Z1, G4double Z2);
    	

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
