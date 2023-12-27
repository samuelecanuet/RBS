class RunManager;
class RunAction;
class DetectorConstruction;
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
#include "G4Material.hh"
#include "G4NistManager.hh"

using namespace std;


class Messenger : public G4UImessenger
{
    
  G4UIcmdWithAString *input_particle;
  G4UIcmdWithADoubleAndUnit *input_energy;
  G4UIcmdWith3VectorAndUnit* input_position;
  G4UIcmdWith3Vector* input_direction;
  G4UIcmdWithAString* input_size;

  G4UIcmdWithAString* input_layer;
  G4UIcmdWithADouble* input_nblayer;
  G4UIcmdWithAString* input_detector;

public:
  Messenger();
  Messenger(RunAction*);
  Messenger(DetectorConstruction*);
  virtual ~Messenger();


  void DefineInputCommands();
  void DefineInputCommandsDetector();
  void DefineInputCommandsRun();

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

  RunAction* fRunAction;
  DetectorConstruction* fDetector;

  vector<pair<G4Material*, G4double>> layers;
  G4int numberOfLayer;
  G4Material* GetMaterial(G4String formula, G4double density);
  std::map<std::string, int> parseChemicalFormula(G4String formula);
  vector<pair<G4Material*, G4double>> GetTarget();

};


