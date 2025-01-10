/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include <vector>
#include <numeric>
#include "TH1D.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction *runAction, Tracking *tracking)
    : G4UserEventAction(),
      fRunAction(runAction), fTracking(tracking)
{   
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *event)
{
    StepVector.clear();

    InitialParticle = event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *event)
{
    
}

void EventAction::AddParticle(G4ThreeVector position, G4ThreeVector direction, G4double energy, G4double time, G4Material *mat)
{
    Data info;
    info.Position = position;
    info.Direction = direction;
    info.KineticEnergy = energy;
    info.Time = time;
    info.Material = mat;

    StepVector.push_back(info);
}

std::vector<Data> EventAction::GetParticleVector()
{
    return StepVector;
}

const G4ParticleDefinition *EventAction::GetInitialParticle()
{
    return InitialParticle;
}
