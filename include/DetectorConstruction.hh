/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"

class G4VPhysicalVolume;
class Messenger;


using namespace std;

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume *Construct();

  G4VPhysicalVolume *Physic_Layer;
  G4VPhysicalVolume *Physics_Detector;
  G4LogicalVolume *logicWorld;
  G4double SuppCatcher_Catcher_radius_inner;

  Messenger *msgg;
  std::vector<std::pair<G4Material *, G4double>> vec;

  void SetTarget(vector<pair<G4Material *, G4double>>);
  void SetDetector(G4String Name, G4Material *Material, G4double Thickness, G4double Area, G4double Angle, G4double Distance, G4double Resolution);

  std::vector<G4Material *> GetMaterials();
  std::vector<G4Material *> vecs;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
