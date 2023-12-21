#include "Messenger.hh"

//----------------------------------------------------------------------
// constructor of the messenger: define the commands

Messenger::Messenger()
{
  DefineInputCommands();
}

// destructor: delete all allocated command objects
Messenger::~Messenger()
{
}

//----------------------------------------------------------------------
// functions defining the commands

// - for input file
void Messenger::DefineInputCommands()
{

  // particle type
  input_particle = new G4UIcmdWithAString("/particle", this);

  // particle energy
  input_energy = new G4UIcmdWithADoubleAndUnit("/energy", this);

  // particle position
  input_position = new G4UIcmdWith3VectorAndUnit("/position", this);

  // particle direction
  input_direction = new G4UIcmdWith3Vector("/direction", this);

  // particle direction
  input_size = new G4UIcmdWithAString("/direction", this);

  
//   // Set Catcher Thickness
//   input_cmd_catcher_thickness = new G4UIcmdWithAString("/Catcher_MylarDeltaThickness", this);
//   input_cmd_catcher_thickness->SetGuidance("Set Catcher thickness");
//   input_cmd_catcher_thickness->SetParameterName("catcher_thickness", false);
//   input_cmd_catcher_thickness->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//----------------------------------------------------------------------
// function processing the commands
void Messenger::SetNewValue(G4UIcommand *cmd, G4String args)
{

  // input file commands

  if (cmd == input_particle)
  {
    particle = G4ParticleTable::GetParticleTable()->FindParticle(args);
  }

  if (cmd == input_energy)
  {
    energy = input_energy->GetNewDoubleValue(args);
  }

  if (cmd == input_position)
  {
    particle_position = input_position->GetNew3VectorValue(args);
  }

  if (cmd == input_direction)
  {
    particle_direction= input_direction->GetNew3VectorValue(args);
  }

  if (cmd == input_size)
  {
    G4double valuex, valuey;
    G4String unit;
    std::istringstream iss(args);
    iss >> valuex >> valuey >> unit;
    particle_size[0] = valuex * G4UnitDefinition::GetValueOf(unit);
    particle_size[1] = valuey * G4UnitDefinition::GetValueOf(unit);
  }
}

G4ParticleDefinition* Messenger::GetParticle(){  return particle;  }

G4double Messenger::GetParticleEnergy(){    return energy;}

G4ThreeVector Messenger::GetParticlePosition(){  return particle_position;}

G4ThreeVector Messenger::GetParticleDirection(){ return particle_direction;}

G4double* Messenger::GetParticleSize(){ return particle_size;}