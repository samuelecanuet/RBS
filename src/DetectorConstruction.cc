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

#include <cmath>

#include "Messenger.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
   
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // World-----------------------------------------------------------------------
  // world parameters

  G4double world_sizeXY = 5*cm;
  G4double world_sizeZ  = 5*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //// CATCHER //////
  /// CENTRAL
  double PEEK_thikness = 1.5 * mm;
  G4double SuppCatcher_Catcher_radius_inner = 7.5 * mm;
  G4double SuppCatcher_Catcher_radius_outer = 9 * mm;
  G4ThreeVector Catcher_Position = G4ThreeVector(0, 0, - PEEK_thikness/2);
  G4Material *Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
  G4Material *Mylar = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");
  G4Material *Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  G4Material *Au = G4NistManager::Instance()->FindOrBuildMaterial("G4_Au");
  thicknessMylarSource = 0.58 * um;
  thicknessAlSource = 50. * nm;

  G4Material *PEEK = new G4Material("PEEKS", 1.32 * g / cm3, 3);
  G4NistManager *man = G4NistManager::Instance();
  PEEK->AddElement(man->FindOrBuildElement("C"), 19);
  PEEK->AddElement(man->FindOrBuildElement("H"), 12);
  PEEK->AddElement(man->FindOrBuildElement("O"), 3);
  G4Tubs *PEEK_Ring = new G4Tubs("PEEK_ring", SuppCatcher_Catcher_radius_inner, SuppCatcher_Catcher_radius_outer, PEEK_thikness / 2, 0., 360 * deg);

  G4LogicalVolume *logicPEEK_Ring = new G4LogicalVolume(PEEK_Ring, PEEK, "logic_Supp_catcher");
  physPEEK_ring = new G4PVPlacement(0, // no rotation
                                            G4ThreeVector(0, 0, + PEEK_thikness / 2),
                                            logicPEEK_Ring, // its fLogical volume
                                            "centralPEEK_Ring",     // its name
                                            logicWorld,   // its mother  volume
                                            false,                  // no boolean operation
                                            0);

  G4VisAttributes *visAtt_PEEK = new G4VisAttributes(G4Colour(1., 1., 1.));
  visAtt_PEEK->SetVisibility(true);
  visAtt_PEEK->SetForceWireframe(false);
  visAtt_PEEK->SetForceSolid(true);
  logicPEEK_Ring->SetVisAttributes(visAtt_PEEK);

  G4Tubs *AlSource1 = new G4Tubs("AlSource1", 0., SuppCatcher_Catcher_radius_inner, thicknessAlSource / 2, 0., 360. * deg);
  G4LogicalVolume *Logic_AlSource1 = new G4LogicalVolume(AlSource1, Al, "LogicAlSource1");               // solid, material, name
  Physics_AlSource1 = new G4PVPlacement(0,                                                                       // no rotation
                                                G4ThreeVector(0., 0., thicknessAlSource / 2), // position
                                                Logic_AlSource1, "LogicAlSource1",                       // its fLogical volume
                                                logicWorld,                                                    // its mother volume
                                                false,                                                                   // no boolean op.
                                                1);                                                                      // copy nb.

  G4Tubs *MylarSource = new G4Tubs("MylarSource", 0., SuppCatcher_Catcher_radius_inner, thicknessMylarSource / 2, 0., 360. * deg);
  G4LogicalVolume *Logic_MylarSource = new G4LogicalVolume(MylarSource, Mylar, "LogicMylarSource");                               // solid, material, name
  Physics_MylarSource = new G4PVPlacement(0,                                                                                                      // no rotation
                                                  G4ThreeVector(0., 0., thicknessMylarSource / 2 + thicknessAlSource), // position
                                                  Logic_MylarSource, "LogicMylarSource",                                                 // its fLogical volume
                                                  logicWorld,                                                                                   // its mother volume
                                                  false,                                                                                                  // no boolean op.
                                                  2);                                                                                                     // copy nb.

  G4Tubs *AlSource2 = new G4Tubs("AlSource2", 0., SuppCatcher_Catcher_radius_inner, thicknessAlSource / 2, 0., 360. * deg);
  G4LogicalVolume *Logic_AlSource2 = new G4LogicalVolume(AlSource2, Al, "LogicAlSource2");                                                                  // solid, material, name
  Physics_AlSource2 = new G4PVPlacement(0,                                                                                                                          // no rotation
                                                G4ThreeVector(0., 0., thicknessAlSource / 2 + thicknessMylarSource + thicknessAlSource), // position
                                                Logic_AlSource2, "LogicAlSource2",                                                                          // its fLogical volume
                                                logicWorld,                                                                                                       // its mother volume
                                                false,                                                                                                                      // no boolean op.
                                                3);                                                                                                                         // copy nb.

  G4VisAttributes *MylarSource_att = new G4VisAttributes(G4Colour(0.94, 0.5, 0.5)); // pink
  MylarSource_att->SetVisibility(true);
  MylarSource_att->SetForceWireframe(false);
  MylarSource_att->SetForceSolid(true);
  Logic_MylarSource->SetVisAttributes(MylarSource_att);

  G4VisAttributes *AlSource_att = new G4VisAttributes(G4Colour(0.94, 0.2, 0.5)); // pink
  AlSource_att->SetVisibility(true);
  AlSource_att->SetForceWireframe(false);
  AlSource_att->SetForceSolid(true);
  Logic_AlSource1->SetVisAttributes(AlSource_att);
  Logic_AlSource2->SetVisAttributes(AlSource_att);

  ///STANDARD CALIB
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

  ///////Detectors
  G4double Active_Surface = 25*mm;
  G4double Camberra_thickness = 100*um;
  G4Tubs *Camberra = new G4Tubs("Camberra", 0., 2.82*mm, Camberra_thickness/2, 0., 360*deg);
  
  ///RBS
  G4RotationMatrix *myRotation = new G4RotationMatrix();
  myRotation->rotateX(45. * deg);
  myRotation->rotateY(0. * deg);
  myRotation->rotateZ(0. * deg);
  G4double RBS_dist = 23*mm;
  G4LogicalVolume *Logic_RBS = new G4LogicalVolume(Camberra, Si, "RBS");                                                                  // solid, material, name
  Physics_RBS = new G4PVPlacement(myRotation,                                                                                                                          // no rotation
                                                G4ThreeVector(0., -RBS_dist*cos(45*deg), -RBS_dist*sin(45*deg)), // position
                                                Logic_RBS, "RBS",                                                                          // its fLogical volume
                                                logicWorld,                                                                                                       // its mother volume
                                                false,                                                                                                                      // no boolean op.
                                                0); 

  ///STIM
  G4double STIM_dist = 5*mm;
  G4LogicalVolume *Logic_STIM = new G4LogicalVolume(Camberra, Si, "STIM");                                                                  // solid, material, name
  Physics_STIM = new G4PVPlacement(0,                                                                                                                          // no rotation
                                                G4ThreeVector(0., 0., STIM_dist), // position
                                                Logic_STIM, "STIM",                                                                          // its fLogical volume
                                                logicWorld,                                                                                                       // its mother volume
                                                false,                                                                                                                      // no boolean op.
                                                0); 
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTXT(G4String name)
{
  TXT = name;
}