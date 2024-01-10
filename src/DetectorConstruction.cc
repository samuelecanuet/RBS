/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolumeStore.hh"

#include <cmath>

#include "Messenger.hh"
#include "Detector.hh"
#include "Layer.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction()
{
  msgg = new Messenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{

  // Get nist material manager
  G4NistManager *nist = G4NistManager::Instance();

  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // World-----------------------------------------------------------------------
  // world parameters

  G4double world_sizeXY = 7 * cm;
  G4double world_sizeZ = 7 * cm;
  G4Material *world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  G4Box *solidWorld =
      new G4Box("World",                                                    // its name
                0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ); // its size

  logicWorld =
      new G4LogicalVolume(solidWorld, // its solid
                          world_mat,  // its material
                          "World");   // its name

  G4VPhysicalVolume *physWorld =
      new G4PVPlacement(0,               // no rotation
                        G4ThreeVector(), // at (0,0,0)
                        logicWorld,      // its logical volume
                        "World",         // its name
                        0,               // its mother  volume
                        false,           // no boolean operation
                        0,               // copy number
                        checkOverlaps);  // overlaps checking

  SuppCatcher_Catcher_radius_inner = 7.5 * mm;

  /// STANDARD CALIB
  // G4double Au_thickness = 10*um;
  // G4Tubs *Au1 = new G4Tubs("Au1", 0., SuppCatcher_Catcher_radius_inner, Au_thickness / 2, 0., 360. * deg);
  // G4LogicalVolume *Logic_Au = new G4LogicalVolume(Au1, Au, "LogicAu");               // solid, material, name
  // G4VPhysicalVolume* phys_Au = new G4PVPlacement(0,                                                                       // no rotation
  //                                               G4ThreeVector(0., 0., Au_thickness / 2), // position
  //                                               Logic_Au, "LogicAu",                       // its fLogical volume
  //                                               logicWorld,                                                    // its mother volume
  //                                               false,                                                                   // no boolean op.
  //                                               0);

  // G4double Si_thickness = 500*um;
  // G4Tubs *Si1 = new G4Tubs("Si1", 0., SuppCatcher_Catcher_radius_inner, Si_thickness / 2, 0., 360. * deg);
  // G4LogicalVolume *Logic_Si = new G4LogicalVolume(Si1, Si, "LogicSi");               // solid, material, name
  // G4VPhysicalVolume *phys_Si = new G4PVPlacement(0,                                                                       // no rotation
  //                                               G4ThreeVector(0., 0., Au_thickness + Si_thickness / 2), // position
  //                                               Logic_Si, "LogicSi",                       // its fLogical volume
  //                                               logicWorld,                                                    // its mother volume
  //                                               false,                                                                   // no boolean op.
  //                                               0);


  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTarget(vector<pair<G4Material *, G4double>> Target)
{

  vector<G4Material *> materialVector;

  G4int Layer_index = 0;
  G4double Position = 0;
  for (const pair<G4Material *, G4double> &layer : Target)
  {
    materialVector.push_back(layer.first);
    Layer_index++;

    std::ostringstream oss;
    oss << "Layer" << Layer_index << "_" << layer.first->GetName();
    G4String name = oss.str();

    G4Tubs *Layer_tubs = new G4Tubs(name, 0., SuppCatcher_Catcher_radius_inner, layer.second / 2, 0., 360. * deg);
    G4LogicalVolume *Logic_Layer = new G4LogicalVolume(Layer_tubs, layer.first, name);
    Physic_Layer = new G4PVPlacement(0, G4ThreeVector(0., 0., Position + layer.second / 2), Logic_Layer, name, logicWorld, false, Layer_index);
    Position += layer.second;

    Logic_Layers.push_back(Logic_Layer);
    Layers.push_back(new Layer(Layer_index - 1, layer.first, fArea, fAngle, fDistance));
    Logic_Layer->SetSensitiveDetector(Layers[Layers.size() - 1]);

    Total_Thickness += layer.second;
  }

  std::set<G4Material *> uniqueMaterials(materialVector.begin(), materialVector.end());
  std::vector<G4Material *> uniqueMaterialVector(uniqueMaterials.begin(), uniqueMaterials.end());
  vecs = uniqueMaterialVector;

  


  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction::SetDetector(G4String Name, G4int Number, G4Material *Material, G4double Thickness, G4double Area, G4double Angle, G4double Distance, G4double Resolution)
{
  G4ThreeVector Position = G4ThreeVector(0., Distance * sin(Angle * deg), Distance * cos(Angle * deg));

  G4RotationMatrix *Rotation = new G4RotationMatrix();
  Rotation->rotateX(Angle * deg);
  Rotation->rotateY(0. * deg);
  Rotation->rotateZ(0. * deg);

  G4Tubs *Detector_tubs = new G4Tubs(Name, 0., sqrt(Area / M_PI), Thickness / 2, 0., 360. * deg);
  G4LogicalVolume *Logic_Detector = new G4LogicalVolume(Detector_tubs, Material, Name);
  Physics_Detector = new G4PVPlacement(Rotation, Position, Logic_Detector, Name, logicWorld, false, 0);

  Distance = Distance - Thickness/2 - 50/2*nm;
  G4Tubs *Detector_tubs1 = new G4Tubs("1", 0., sqrt(Area / M_PI), 50*nm/2, 0., 360. * deg);
  G4LogicalVolume *Logic_Detector1 = new G4LogicalVolume(Detector_tubs1, Material, "1");
  G4VPhysicalVolume* Physics_Detector1 = new G4PVPlacement(Rotation, G4ThreeVector(0., Distance * sin(Angle * deg), Distance * cos(Angle * deg)), Logic_Detector1, "1", logicWorld, false, 0);

  Detectors.push_back(new Detector(Name, Number, Resolution));
  Logic_Detector->SetSensitiveDetector(Detectors[Detectors.size() - 1]);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  if (Name == "RBS")
  {
    fArea = Area;
    fDistance = Distance;
    fAngle = Angle;
  }
}

std::vector<G4Material *> DetectorConstruction::GetMaterials()
{
  return vecs;
}

std::vector<Detector *> DetectorConstruction::GetDetectors()
{
  return Detectors;
}

std::vector<Layer *> DetectorConstruction::GetLayers()
{
  return Layers;
}

G4double DetectorConstruction::GetTotalThickness()
{
  return Total_Thickness;
}