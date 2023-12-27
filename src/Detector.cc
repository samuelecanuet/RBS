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
//----------------------------------------------------------------------

Detector::Detector()
    : G4VSensitiveDetector("")
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

  G4cout<<"HIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIT"<<G4endl;

  return (true);
}
