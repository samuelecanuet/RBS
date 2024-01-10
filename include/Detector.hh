
#include "G4VSensitiveDetector.hh"
#include "G4EventManager.hh"
#include "globals.hh"


class RunAction;
using namespace CLHEP;

class Detector : public G4VSensitiveDetector
{

  //------------------------------------------------------------
public:
  Detector(G4String, G4int, G4double);
  ~Detector();

  void Initialize(G4HCofThisEvent *);
  G4bool ProcessHits(G4Step *, G4TouchableHistory *);

  G4EventManager* fEventManager;

  G4String fName;
  G4String GetName();
  G4int fNumber; 
  G4double fResolution;

  G4double Ekin_res;
  RunAction* fRunAction;
};
