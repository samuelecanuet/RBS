/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4Step.hh"

class RunAction;
class Tracking;

struct Data
{
  G4int TrackID;
  G4ThreeVector Position;
  G4ThreeVector Direction;
  G4double KineticEnergy;
  G4double Time;
  G4Material *Material;
};

class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction *runAction, Tracking *tracking);
  virtual ~EventAction();

  virtual void BeginOfEventAction(const G4Event *event);
  virtual void EndOfEventAction(const G4Event *event);
  std::vector<Data> StepVector;

  void AddParticle(G4ThreeVector position, G4ThreeVector direction, G4double energy, G4double time, G4Material *mat);
  std::vector<Data> GetParticleVector();
  const G4ParticleDefinition *InitialParticle;
  const G4ParticleDefinition *GetInitialParticle();

private:
  RunAction *fRunAction;
  Tracking *fTracking;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
