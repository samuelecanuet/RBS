
#include "G4VSensitiveDetector.hh"
#include "G4EventManager.hh"
#include "globals.hh"
#include "G4EmCalculator.hh"

class DetectorConstruction;
class Messenger;
class RunAction;

using namespace CLHEP;


class Layer : public G4VSensitiveDetector
{

  //------------------------------------------------------------
public:
  Layer(G4int, G4Material *, G4double, G4double, G4double);
  ~Layer();

  void Initialize(G4HCofThisEvent *);

  G4bool ProcessHits(G4Step *, G4TouchableHistory *);

  G4EventManager *fEventManager;
  RunAction* fRunAction;

  G4int GetNumber();
  G4int fNumber;

  DetectorConstruction *fdet;
  G4double CollisionEnergy(G4double Energy, const G4Element *Element, G4double Angle);
  G4ThreeVector DirectionRBS();
  G4double CrossSectionRBS(G4double Energy, const G4Element *Element, G4double theta);
  const G4Element *RandomElement(G4Material *Material);


  G4double GetRandomTheta(G4double a);

  G4double CalcAngleCMFrame(G4double angle, G4double M1, G4double M2);
  G4double CalcEnergyCMFrame(G4double energy, G4double M1, G4double M2);
  G4double CalcDiffRuthXsecLAB(G4double M1, G4double M2, G4double angle, G4double xsection);
  G4double CalcDiffRuthXsecCM(G4double E, G4double angleCM, G4double Z1, G4double Z2);
  G4double CalculateTotalRBSYield(G4double energy, G4double M1, G4double M2, G4double Z1, G4double Z2, G4double angle);
  G4double CalcAndersenScreening(G4double energy_cm, G4double angle_cm, G4double Z1, G4double Z2);

  G4double CalcNuclEnStraggling(G4double Z1, G4double Z2, G4double M1, G4double M2, G4double atdens,G4double distance);
  G4double CalcBohrStrag(G4double Z1, G4double Z2, G4double atomDens, G4double dist);
  G4double CalculateTotalBohrStraggling(G4double energy, G4ParticleDefinition* particle, G4Material* mat, G4double distance);


  const G4ParticleDefinition *InitialParticle;
  std::vector<G4double> energy, step_lenght;
  std::vector<G4ThreeVector> position;
  std::vector<G4ThreeVector> direction;


  G4EmCalculator emCalculator;
  G4Material *fMaterial;
  G4double fArea, fAngle, fDistance;

  Messenger* msg;
};
