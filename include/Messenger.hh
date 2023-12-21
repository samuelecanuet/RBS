#ifndef WITCH_MESSENGER_HH
#define WITCH_MESSENGER_HH


class RunManager; // just indicates that this class exists

// Geant4 include files
#include "G4UImessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UnitsTable.hh"

#include "G4ParticleTable.hh"


class Messenger : public G4UImessenger
{
    
  G4UIcmdWithAString *input_particle;
  G4UIcmdWithADoubleAndUnit *input_energy;
  G4UIcmdWith3VectorAndUnit* input_position;
  G4UIcmdWith3Vector* input_direction;
  G4UIcmdWithAString* input_size;

public:
  Messenger();
  virtual ~Messenger();


  void DefineInputCommands();
  void SetNewValue(G4UIcommand *cmd, G4String args);

  G4ParticleDefinition* particle;
  G4ParticleDefinition* GetParticle();

  G4double energy;
  G4double GetParticleEnergy();

  G4ThreeVector particle_position;
  G4ThreeVector GetParticlePosition();

  G4ThreeVector particle_direction;
  G4ThreeVector GetParticleDirection();

  G4double particle_size[2];
  G4double* GetParticleSize();


};

#endif
