#include "DetectorConstruction.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();

  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  G4bool isotopes = false;
  G4String name, symbol;
  G4int ncomponents, natoms;
  G4double fractionmass;
  
  G4Element* H = nistManager->FindOrBuildElement("H", isotopes);
  G4Element* N = nistManager->FindOrBuildElement("N", isotopes);
  G4Element* C = nistManager->FindOrBuildElement("C", isotopes);
  G4Element* O = nistManager->FindOrBuildElement("O", isotopes);
  G4Element* Ar = nistManager->FindOrBuildElement("Ar", isotopes);
  
  G4Element* W = new G4Element("Wolfram","W", 74., 183.84*g/mole);
  G4Element* Cr = new G4Element("Chromium","Cr", 24., 52.00*g/mole);
  G4Element* Fe = new G4Element("Ferrum","Fe", 26., 55.85*g/mole);
  G4Element* Ni = new G4Element("Nickel","Ni", 28., 58.70*g/mole);
  G4Element* Ti = new G4Element("Titanium", "Ti", 22., 47.88*g/mole);
  G4Element* Xe = new G4Element("Xenon", "Xe", 54., 131.29*g/mole);

  G4Material* WMat =
      new G4Material("WMat", density = 19.25 * CLHEP::g / CLHEP::cm3,
                     ncomponents = 1, kStateSolid, 293.*kelvin, 2.*atmosphere);
                     
  WMat->AddElement(W, fractionmass = 1.0);

  G4Material* ArXe_70_30 =
      new G4Material("ArXe_70_30", density = 3.0142 * CLHEP::g / CLHEP::cm3,
                     ncomponents = 2, kStateGas, 293.*kelvin, 2.*atmosphere);
                     
  ArXe_70_30->AddElement(Ar, fractionmass = 0.7);
  ArXe_70_30->AddElement(Xe, fractionmass = 0.3);

  auto steel_mat = 
  new G4Material("steel_mat", 7.9 * g / cm3, ncomponents=5, kStateSolid, 293.*kelvin, 2.*atmosphere);
  
    steel_mat->AddElement(C, fractionmass=0.008);
    steel_mat->AddElement(Cr, fractionmass=0.18);
  
    steel_mat->AddElement(Fe, fractionmass=0.707);
    steel_mat->AddElement(Ni, fractionmass=0.10);
    steel_mat->AddElement(Ti, fractionmass=0.005);
  
  //Air material
  G4Material* worldmaterial =
      new G4Material("worldmaterial", density = 1.290 * CLHEP::g / CLHEP::cm3,
                     ncomponents = 2, kStateGas, 293.*kelvin, 1.*atmosphere);
                     
  worldmaterial->AddElement(N, fractionmass = 0.7);
  worldmaterial->AddElement(O, fractionmass = 0.3);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
 // G4double absoThickness = 1.0 * m;
 // G4double gapThickness = 5.0 * mm;

  auto worldSizeXYZ = 1.0 * m;

  G4double wireRadius = 0.05 * mm;
  G4double tubeRadius = 29 * mm;
  G4double tubeHalfLength = 150 * mm;
  G4double tubeThickness = 0.5 * mm;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("worldmaterial");
  auto absorberMaterial = G4Material::GetMaterial("worldmaterial");
  auto gapMaterial = G4Material::GetMaterial("ArXe_70_30");

  //auto steelMaterial = G4Material::GetMaterial("steel_mat");

   auto wireMaterial = G4Material::GetMaterial("WMat");

  if (!defaultMaterial || !absorberMaterial || !gapMaterial || !wireMaterial) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()", "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS = new G4Box("World",  // its name
                          worldSizeXYZ / 2, worldSizeXYZ / 2, worldSizeXYZ / 2);  // its size

  auto worldLV = new G4LogicalVolume(worldS,  // its solid
                                     defaultMaterial,  // its material
                                     "World");  // its name

  auto worldPV = new G4PVPlacement(nullptr,  // no rotation
                                   G4ThreeVector(),  // at (0,0,0)
                                   worldLV,  // its logical volume
                                   "World",  // its name
                                   nullptr,  // its mother  volume
                                   false,  // no boolean operation
                                   0,  // copy number
                                   fCheckOverlaps);  // checking overlaps
  //
  // Absorber
  //
  auto absorberS = new G4Box("Abso",  // its name
                             worldSizeXYZ / 3, worldSizeXYZ / 3, worldSizeXYZ / 3);  // its size

  auto absorberLV = new G4LogicalVolume(absorberS,  // its solid
                                        defaultMaterial,  // its material
                                        "Abso");  // its name

  fAbsorberPV = new G4PVPlacement(nullptr,  // no rotation
                                  G4ThreeVector(0., 0., 0.),  // its position
                                  absorberLV,  // its logical volume
                                  "Abso",  // its name
                                  worldLV,  // its mother  volume
                                  false,  // no boolean operation
                                  0,  // copy number
                                  fCheckOverlaps);  // checking overlaps

    // Wire
  //
  G4VSolid* wireS =
      new G4Tubs("Wire",                                     // its name
                 0, wireRadius, tubeHalfLength, 0, 2 * pi);  // its size

  G4LogicalVolume* wireLV = new G4LogicalVolume(wireS,         // its solid
                                                gapMaterial,  // its material
                                                "Wire");       // its name

  G4VPhysicalVolume* WirePV = new G4PVPlacement(0,                          // no rotation
                              G4ThreeVector(),  // its position
                              wireLV,                     // its logical volume
                              "Wire",                     // its name
                              worldLV,                     // its mother  volume
                              false,            // no boolean operation
                              0,                // copy number
                              fCheckOverlaps);  // checking overlaps
  
  //
  // Gap
  //
  auto gapS = new G4Tubs("Gas",  // its name
                              wireRadius, tubeRadius,
                              tubeHalfLength, 0, 2 * pi);
                              
  auto gapLV = new G4LogicalVolume(gapS,  // its solid
                                   gapMaterial,  // its material
                                   "Gap");  // its name

  fGapPV = new G4PVPlacement(0,  // no rotation
                             G4ThreeVector(),  // its position
                             gapLV,  // its logical volume
                             "Gap",  // its name
                             worldLV,  // its mother  volume
                             false,  // no boolean operation
                             0,  // copy number
                             fCheckOverlaps);  // checking overlaps

  // Outer tube
  //
  /*G4VSolid* tubeS = new G4Tubs("Tube",  // its name
                               tubeRadius, tubeRadius + tubeThickness, tubeHalfLength + tubeThickness, 0,
                               2 * pi);  // its size

  G4LogicalVolume* tubeLV =
      new G4LogicalVolume(tubeS,            // its solid
                          steelMaterial,  // its material
                          "Tube");          // its name


  G4VPhysicalVolume* TubePV = new G4PVPlacement(
      0,
      G4ThreeVector(),  // its position
      tubeLV,                                          // its logical volume
      "Tube",                                          // its name
      worldLV,                                         // its mother  volume
      false,                                           // no boolean operation
      0,                                               // copy number
      fCheckOverlaps);*/                                 // checking overlaps
      
  
  // Visualization attributes

  G4VisAttributes* VisAttYellow = new G4VisAttributes(G4Colour::Yellow());
  G4VisAttributes* VisAttBlue = new G4VisAttributes(G4Colour::Blue());
  G4VisAttributes* VisAttGreen = new G4VisAttributes(G4Colour::Green());
  G4VisAttributes* VisAttRed = new G4VisAttributes(G4Colour::Red());
  G4VisAttributes* VisAttWhite = new G4VisAttributes(G4Colour::White());
  
  VisAttWhite->SetVisibility(true);
  VisAttWhite->SetForceWireframe(true);
  worldLV->SetVisAttributes(VisAttWhite);
  
  VisAttRed->SetVisibility(true);
  VisAttRed->SetForceWireframe(true);
  wireLV->SetVisAttributes(VisAttRed);
  
 /*VisAttGreen->SetVisibility(true);
  VisAttGreen->SetForceWireframe(true);
  tubeLV->SetVisAttributes(VisAttGreen);*/
  
  //gapLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  VisAttYellow->SetVisibility(true);
  VisAttYellow->SetForceWireframe(true);
  gapLV->SetVisAttributes(VisAttYellow);


  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B4
