#include "Detector.hh"
#include "G4VProcess.hh"
#include "G4TrackingManager.hh"
#include "G4Event.hh"
#include <G4Field.hh>
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4LogicalVolume.hh"
#include "G4TouchableHandle.hh"
#include "G4Navigator.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "globals.hh"

//----------------------------------------------------------------------

Detector::Detector(G4String Name, G4int Number, G4double Resolution)
    : G4VSensitiveDetector(Name), fName(Name), fNumber(Number), fResolution(Resolution)
{
}

// destructor
Detector::~Detector()
{
}

//----------------------------------------------------------------------

void Detector::Initialize(G4HCofThisEvent *)
{
}

G4bool Detector::ProcessHits(G4Step *step, G4TouchableHistory *)
{
  if ( step->IsFirstStepInVolume() && G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition() == step->GetTrack()->GetParticleDefinition() )
  {

    G4double Ekin = step->GetPreStepPoint()->GetKineticEnergy() / keV;
    Ekin_res = Ekin + G4RandGauss::shoot(0, fResolution/keV);

    fRunAction->FillNtupleDColumn(1, fNumber*3, Ekin);
    fRunAction->FillNtupleDColumn(1, fNumber*3+1, Ekin_res);
    fRunAction->FillNtupleDColumn(1, fNumber*3+2, step->GetTrack()->GetWeight());
    fRunAction->AddNtupleRow(1);

    fRunAction->FillH1(3*fNumber, Ekin_res);
    fRunAction->FillH1(3*fNumber+1, Ekin_res,  step->GetTrack()->GetWeight());
    fRunAction->FillH1(3*fNumber+2, Ekin_res,  1/step->GetTrack()->GetWeight());

    step->GetTrack()->SetTrackStatus(fStopAndKill);
  }

  return (true);
}

G4String Detector::GetName()
{
  return fName;
}