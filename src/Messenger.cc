#include "Messenger.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"

//----------------------------------------------------------------------
// constructor of the messenger: define the commands

Messenger::Messenger()
{
  DefineInputCommands();
}

Messenger::Messenger(RunAction *runAction) : fRunAction(runAction)
{
  DefineInputCommandsRun();
}

Messenger::Messenger(DetectorConstruction *detectorConstruction) : fDetector(detectorConstruction)
{
  DefineInputCommandsDetector();
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

  // particle size
  input_size = new G4UIcmdWithAString("/size", this);
}

void Messenger::DefineInputCommandsDetector()
{
  // nblayer
  input_nblayer = new G4UIcmdWithADouble("/numberOfLayer", this);

  // layer
  input_layer = new G4UIcmdWithAString("/layer", this);

   // detector
  input_detector = new G4UIcmdWithAString("/detector", this);
}

void Messenger::DefineInputCommandsRun()
{
   // particle energy
  input_energy = new G4UIcmdWithADoubleAndUnit("/energy", this);
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
    particle_direction = input_direction->GetNew3VectorValue(args);
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

   if (cmd == input_layer)
  {
    G4int number;
    G4double thickness;
    G4String material, unit;

    std::istringstream iss(args);
    iss >> number >> material >> thickness >> unit;

    layers.push_back(make_pair(G4NistManager::Instance()->FindOrBuildMaterial(material), thickness * G4UnitDefinition::GetValueOf(unit)));

    if (layers.size() == numberOfLayer)
    {
      fDetector->SetTarget(layers);
    }
  }

  if (cmd == input_nblayer)
  {
    numberOfLayer = input_nblayer->GetNewDoubleValue(args);
  }

  if (cmd == input_detector)
  {
    G4double thickness, distance, angle, resolution, area;
    G4String material, thickness_unit, distance_unit, name, resolution_unit, area_unit;

    std::istringstream iss(args);
    iss >> name >> material >> thickness >> thickness_unit >> area >> area_unit >> angle >> distance >> distance_unit >> resolution >> resolution_unit;

    fDetector->SetDetector(name, G4NistManager::Instance()->FindOrBuildMaterial(material), thickness * G4UnitDefinition::GetValueOf(thickness_unit), area * G4UnitDefinition::GetValueOf(area_unit), angle, distance * G4UnitDefinition::GetValueOf(distance_unit), resolution * G4UnitDefinition::GetValueOf(resolution_unit));
  }

}

G4ParticleDefinition *Messenger::GetParticle() { return particle; }
G4double Messenger::GetParticleEnergy() { return energy; }
G4ThreeVector Messenger::GetParticlePosition() { return particle_position; }
G4ThreeVector Messenger::GetParticleDirection() { return particle_direction; }
G4double *Messenger::GetParticleSize() { return particle_size; }
vector<pair<G4Material*, G4double>> Messenger::GetTarget(){return layers;}
