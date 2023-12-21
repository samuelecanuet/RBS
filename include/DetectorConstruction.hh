/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"

class G4VPhysicalVolume;
class Messenger;


/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

  G4VPhysicalVolume *Physics_AlSource1;
  G4VPhysicalVolume *Physics_MylarSource;
  G4VPhysicalVolume *Physics_AlSource2;
  G4VPhysicalVolume *physPEEK_ring;
  G4VPhysicalVolume *Physics_RBS;
  G4VPhysicalVolume *Physics_STIM;

  G4double thicknessMylarSource;
  G4double thicknessAlSource;

  Messenger* msg;

  void SetTXT(G4String name);
  G4String TXT = "det non";


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

