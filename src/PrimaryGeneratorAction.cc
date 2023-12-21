/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "Messenger.hh"
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
{
  fMessenger = new Messenger(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{

  
  gun.SetParticleDefinition(        fMessenger->GetParticle());
  gun.SetParticleEnergy(            fMessenger->GetParticleEnergy());
  gun.SetParticlePosition(          RandomXY(fMessenger->GetParticlePosition(), fMessenger->GetParticleSize()));
  gun.SetParticleMomentumDirection( fMessenger->GetParticleDirection());
  gun.GeneratePrimaryVertex(        anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition* PrimaryGeneratorAction::GetInitialParticle()
{
  return fMessenger->GetParticle();
}

G4ThreeVector PrimaryGeneratorAction::RandomXY(G4ThreeVector Position, G4double xy[])
{
  Position.x += G4RandGauss::shoot(Position.x, xy[0]);
  Position.y += G4RandGauss::shoot(Position.y, xy[1]);
  return Position;
}
