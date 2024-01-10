
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "Messenger.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "Tracking.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include <cmath>
#include "G4CrossSectionDataStore.hh"
#include "G4NistManager.hh"
#include "G4HadronicProcessStore.hh"

#include "G4VProcess.hh"
#include "G4NistManager.hh"
#include "G4ProcessVector.hh"
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction *eventAction, RunAction *runAction, Tracking *tracking, PrimaryGeneratorAction *generator)
	: G4UserSteppingAction(), fEventAction(eventAction), fRunAction(runAction), fTracking(tracking), fGenerator(generator)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step *step)
{

	if (step->GetTrack()->GetDefinition()->GetParticleName() == "e-")
	{
		step->GetTrack()->SetTrackStatus(fStopAndKill);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

