//My
#include "DetectorConstruction.hh"
#include "SensitiveDetector.hh"

//G4
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Paraboloid.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Color.hh"
#include "G4TwoVector.hh"
#include "G4SDManager.hh"
#include "globals.hh"
//magnetic field
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4UserLimits.hh"
//GDML
//#include <G4GDMLParser.hh>

//root 
#include "TMath.h"

DetectorConstruction::DetectorConstruction()
{
  //magField = new MagneticField();
  worldVisAtt = new G4VisAttributes();
  quartzVisAtt = new G4VisAttributes();
  sensitiveVisAtt = new G4VisAttributes();
  pmtboxVisAtt = new G4VisAttributes();
  absVisAtt = new G4VisAttributes();
  // Define Materials to be used
  DefineMaterials();
}

DetectorConstruction::~DetectorConstruction()
{
  //delete magField;
  delete worldVisAtt;
  delete quartzVisAtt;
  delete sensitiveVisAtt;
  delete pmtboxVisAtt;
  delete absVisAtt;
  delete stepLimit;
}

void DetectorConstruction::DefineMaterials()
{
  G4String symbol;
  G4double a, z, density;
  G4int ncomponents, natoms;
  G4double fractionmass;

  // Define elements
  //G4Element* H = 
  //new G4Element("Hydrogen", symbol = "H", z = 1., a = 1.01*g/mole);
  G4Element* C = 
    new G4Element("Carbon",   symbol = "C", z = 6., a = 12.01*g/mole);
  G4Element* N = 
    new G4Element("Nitrogen", symbol = "N", z = 7., a = 14.01*g/mole);
  G4Element* O =
    new G4Element("Oxygen",   symbol = "O", z = 8., a = 16.00*g/mole);
  G4Element* Si = 
    new G4Element("Silicon",  symbol = "Si", z = 14., a = 28.09*g/mole);
  G4Element* Al = 
    new G4Element("Aluminum", symbol = "Al", z = 13., a = 26.98*g/mole);

  // Quartz Material (SiO2)
  G4Material* SiO2 = 
    new G4Material("quartz", density = 2.200*g/cm3, ncomponents = 2);
  SiO2->AddElement(Si, natoms = 1);
  SiO2->AddElement(O , natoms = 2);

  // Quartz Material (SiO2_cladd)
  G4Material* SiO2_cladd = 
    new G4Material("quartzCladd", density = 2.200*g/cm3, ncomponents = 2);
  SiO2_cladd->AddElement(Si, natoms = 1);
  SiO2_cladd->AddElement(O , natoms = 2);

  // Quartz Material (SiO2_coat)
  G4Material* SiO2_coat = 
    new G4Material("quartzCoat", density = 2.200*g/cm3, ncomponents = 2);
  SiO2_coat->AddElement(Si, natoms = 1);
  SiO2_coat->AddElement(O , natoms = 2);
  
  // Air
  G4Material* Air = 
    new G4Material("Air", density = 1.290*mg/cm3, ncomponents = 2);
  //G4Material* Air = 
  //new G4Material("Air", density = 0.000290*mg/cm3, ncomponents = 2);
  Air->AddElement(N, fractionmass = 0.7);
  Air->AddElement(O, fractionmass = 0.3);

  // Aluminum
  G4Material* Aluminum =
    new G4Material("Aluminum", density = 2.7*g/cm3, ncomponents = 1);
  Aluminum->AddElement(Al, fractionmass = 1.0);

  // Assign Materials
  world.material = Air;
  sec1.material = SiO2;
  sec3.material = SiO2;
  fiberCorr.material = SiO2;
  fiberClad.material = SiO2_cladd;
  //fiberCoat.material = SiO2_coat;
  fiberCoat.material = Aluminum;
  sensitive.material = Aluminum;
  abs1.material = Aluminum;
  abs21.material = Aluminum;
  abs22.material = Aluminum;
  abs23.material = Aluminum;
  abs24.material = Aluminum;

  //
  // Generate and Add Material Properties Table
  //						
  const G4int num = 36;
  G4double WaveLength[num];
  G4double Absorption[num];      // Default value for absorption
  G4double AirAbsorption[num];
  G4double AirRefractiveIndex[num];
  G4double PhotonEnergy[num];

  // Absorption of quartz per 1m
  G4double QuartzAbsorption[num] =
    {0.999572036,0.999544661,0.999515062,0.999483019,0.999448285,
     0.999410586,0.999369611,0.999325013,0.999276402,0.999223336,
     0.999165317,0.999101778,0.999032079,0.998955488,0.998871172,
     0.998778177,0.99867541 ,0.998561611,0.998435332,0.998294892,
     0.998138345,0.997963425,0.997767484,0.997547418,0.99729958 ,
     0.99701966 ,0.99670255 ,0.996342167,0.995931242,0.995461041,
     0.994921022,0.994298396,0.993577567,0.992739402,0.991760297,
     0.990610945};

  for (int i=0; i<num; i++) {
    WaveLength[i] = (300 + i*10)*nanometer;
    Absorption[i] = 100*m;      // Fake number for no absorption
    //initially
    AirAbsorption[i] = 4.*cm;   // If photon hits air, kill it
    //AirAbsorption[i] = 1.*cm;   // If photon hits air, kill it
    //AirAbsorption[i] = 0.1*cm;   // If photon hits air, kill it
    //AirAbsorption[i] = 100.0*m;
    AirRefractiveIndex[i] = 1.;
    PhotonEnergy[num - (i+1)] = twopi*hbarc/WaveLength[i];
    /* Absorption is given per length and G4 needs mean free path
       length, calculate it here
       mean free path length - taken as probablility equal 1/e
       that the photon will be absorbed */
    QuartzAbsorption[i] = (-1)/log(QuartzAbsorption[i])*100*cm;
    //EpotekAbsorption[i] = (-1)/log(EpotekAbsorption[i])*
    //epotekBarJoint.thickness;
  }

  G4double QuartzRefractiveIndex[num] =
    {1.456535,1.456812,1.4571  ,1.457399,1.457712,1.458038,
     1.458378,1.458735,1.459108,1.4595  ,1.459911,1.460344,
     1.460799,1.46128 ,1.461789,1.462326,1.462897,1.463502,
     1.464146,1.464833,1.465566,1.46635 ,1.46719 ,1.468094,
     1.469066,1.470116,1.471252,1.472485,1.473826,1.475289,
     1.476891,1.478651,1.480592,1.482739,1.485127,1.487793};

  G4double CladdingRefractiveIndex[num];

  for(int i=0; i<num; i++){
    CladdingRefractiveIndex[i] = TMath::Sqrt(QuartzRefractiveIndex[i]*QuartzRefractiveIndex[i]-0.22*0.22); 
    //CladdingRefractiveIndex[i] = 1.0;
  }

  // Assign absorption and refraction to materials

  // Quartz
  G4MaterialPropertiesTable* QuartzMPT = new G4MaterialPropertiesTable();
  QuartzMPT->AddProperty("RINDEX", PhotonEnergy, QuartzRefractiveIndex, num);
  QuartzMPT->AddProperty("ABSLENGTH", PhotonEnergy, QuartzAbsorption, num);
  
  // Cladding (of the fiber) only for the fiber aplication
  G4MaterialPropertiesTable* CladdingMPT = new G4MaterialPropertiesTable();
  CladdingMPT->AddProperty("RINDEX", PhotonEnergy, CladdingRefractiveIndex, num);
  CladdingMPT->AddProperty("ABSLENGTH", PhotonEnergy, QuartzAbsorption, num);

  // Assign this material to the bars
  sec1.material->SetMaterialPropertiesTable(QuartzMPT);
  fiberCorr.material->SetMaterialPropertiesTable(QuartzMPT);
  fiberClad.material->SetMaterialPropertiesTable(CladdingMPT);

  // Air
  G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
  AirMPT->AddProperty("RINDEX", PhotonEnergy, AirRefractiveIndex, num);
  AirMPT->AddProperty("ABSLENGTH", PhotonEnergy, AirAbsorption, num);
  
  // Assign these properties to the world volume
  world.material->SetMaterialPropertiesTable(AirMPT);
}

G4VPhysicalVolume* DetectorConstruction::Construct(){

  G4double world_sizeX = 40.0*cm;
  G4double world_sizeY = 20.0*cm;
  //Without fibers
  //G4double world_sizeZ = 40.0*cm;
  //With fibers
  G4double world_sizeZ = 80.0*cm;

  G4double sec1_sizeX = 0.5*cm;
  G4double sec1_sizeY = 1.0*cm;
  G4double sec1_sizeZ = 36.63*cm;

  //Volume for subtraction
  G4double sec5_sizeX = sec1_sizeX/TMath::Sin(UA9Const::angleDet);
  G4double sec5_sizeY = sec1_sizeY*2.0;
  G4double sec5_sizeZ = sec1_sizeX;

  bool sec4volIsOn = false;
  G4double sec4_sizeX = 3.0*cm;
  G4double sec4_sizeY = sec1_sizeY;  
  G4double sec4_sizeZ = 5.0*mm;

  bool abs2volIsOn = false;
  G4double abs2_sizeX = 1.0*mm;
  G4double abs2_sizeY = sec1_sizeY;
  G4double abs2_sizeZ = sec4_sizeZ;

  G4double abs3_sizeX = 1.0*mm;
  G4double abs3_sizeY = sec1_sizeY;  
  G4double abs3_sizeZ = 2.0*cm;

  G4double abs4_sizeX = abs3_sizeX;
  G4double abs4_sizeY = abs2_sizeZ;  
  G4double abs4_sizeZ = abs3_sizeZ;

  //G4double xCopyShift = (sec1_sizeX + 0.5)/TMath::Cos(UA9Const::angleDet);
  //Last development
  G4double distBetweenBars = 1.0*mm;
  G4double xCopyShift = -sec1_sizeX - distBetweenBars;
  G4double zCopyShift =  sec1_sizeX/TMath::Tan(UA9Const::angleDet) + distBetweenBars/TMath::Tan(UA9Const::angleDet);

  G4double  para_alpha = 0.0;
  G4double  para_theta = UA9Const::angleDet;
  G4double  para_phi = 0.0;
  G4double  para_dx = sec1_sizeX/TMath::Cos(para_theta)/2.0;
  G4double  para_dy = sec1_sizeY/2.0;
  G4double  para_dz = sec1_sizeZ*TMath::Cos(para_theta)/2.0;

  //in case of 90 deg angle 
  //G4double senDet_sizeX = sec1_sizeX/2.0;
  //G4double senDet_sizeY = sec1_sizeY/2.0;
  //in case of inclination angle 
  //G4double senDet_sizeX = para_dx;
  //G4double senDet_sizeY = para_dy;
  //G4double senDet_sizeZ = 1.0*mm;
  G4double senDet_sizeX = sec5_sizeX;
  G4double senDet_sizeY = sec1_sizeY;
  G4double senDet_sizeZ = 1.0*mm;

  //in case of 90 deg angle 
  //G4double abs1_sizeX = sec1_sizeX/2.0;
  //G4double abs1_sizeY = sec1_sizeY/2.0;
  //in case of inclination angle 
  G4double abs1_sizeX = para_dx;
  G4double abs1_sizeY = para_dy;
  G4double abs1_sizeZ = 1.0*mm;

  G4double fiberCorr_Rmin = 0.0;
  G4double fiberCorr_Rmax = 0.6*mm/2.0;
  G4double fiberCorr_L    = 50*cm;
  G4double fiber_sizeXshift = 0.0*mm;

  G4double fiberClad_Rmin = fiberCorr_Rmax;
  G4double fiberClad_Rmax = 0.66*mm/2.0;
  G4double fiberClad_L    = fiberCorr_L;

  G4double fiberCoat_Rmin = fiberClad_Rmax;
  G4double fiberCoat_Rmax = 0.70*mm/2.0;
  G4double fiberCoat_L    = fiberCorr_L;

  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;

  // 
  // Define World Volume
  //
  world.solid = new G4Box("World",world_sizeX/2.0,world_sizeY/2.0,world_sizeZ/2.0);
  world.logical = new G4LogicalVolume(world.solid,world.material,"World");
  world.physical = new G4PVPlacement(0,G4ThreeVector(),world.logical,"World",0,false,0);

  //
  // Sensitive volume
  //
  sensitive.solid = new G4Box("Sensitive", senDet_sizeX/2.0, senDet_sizeY/2.0, senDet_sizeZ/2.0);
  sensitive.logical = new G4LogicalVolume(sensitive.solid, sensitive.material,"Sensitive");
  //in case of the inclination angle
  //Ta.setX(sec1_sizeZ/2.0*TMath::Sin(para_theta));
  // in case of the 90 deg angle
  //Ta.setX(0.0);
  //Ta.setY(0.0);
  // in case of the 90 deg angle
  //Ta.setZ(sec1_sizeZ/2.0 + senDet_sizeZ/2.0);
  //in case of the inclination angle and fibers
  //Ta.setZ(sec1_sizeZ/2.0*TMath::Cos(para_theta) + fiberCorr_L + senDet_sizeZ/2.0);
  //in case of the inclination angle and no fibers
  //Ta.setZ(sec1_sizeZ/2.0*TMath::Cos(para_theta) + senDet_sizeZ/2.0);
  /////////////////////////
  G4cout<<"MyCalculation"<<G4endl;
  G4RotationMatrix RaMysen;
  G4ThreeVector TaMysen;
  G4Transform3D TrMysen;
  G4double alphaRotsen = twopi/4.0 - UA9Const::angleDet;
  TaMysen.set(senDet_sizeX/2.0,0.0,senDet_sizeZ/2.0);
  //G4double zOld = TaMy.getZ() - sec1_sizeZ/2.0;
  //G4double xOld = TaMy.getX();
  G4cout<<" x "<<TaMysen.getX()<<G4endl
	<<" y "<<TaMysen.getY()<<G4endl
	<<" z "<<TaMysen.getZ()<<G4endl;
  TaMysen.rotateY(alphaRotsen);
  G4cout<<" x "<<TaMysen.getX()<<G4endl
	<<" y "<<TaMysen.getY()<<G4endl
	<<" z "<<TaMysen.getZ()<<G4endl;
  G4double zNewsen = TaMysen.getZ();
  G4double xNewsen = TaMysen.getX();
  /////////////////////////
  Ta.setX( sec1_sizeX/2.0 - xNewsen);
  Ta.setY(0.0);
  Ta.setZ(-sec1_sizeZ/2.0 - zNewsen);
  Ra.rotateY(alphaRotsen);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *sensitivePhysical_zp = new G4PVPlacement(Tr,                //Transformation
  							      sensitive.logical, //its logical volume	 
  							      "Sensitive",       //its name
  							      world.logical,     //its mother volume
  							      false,	         //no boolean operation
  							      0);	         //copy number

  //1
  Ta.setX( sec1_sizeX/2.0 - xNewsen + xCopyShift);
  Ta.setY(0.0);
  Ta.setZ(-sec1_sizeZ/2.0 - zNewsen + zCopyShift);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *sensitivePhysical_zp1 = new G4PVPlacement(Tr,               //Transformation
  							      sensitive.logical, //its logical volume	 
  							      "Sensitive",       //its name
  							      world.logical,     //its mother  volume
  							      false,	         //no boolean operation
  							      0);	         //copy number
  Ra.rotateY(-alphaRotsen);
  //2
  Ta.setX(sec1_sizeZ/2.0*TMath::Sin(para_theta) + 2*xCopyShift);
  Ta.setY(0.0);
  Ta.setZ(sec1_sizeZ/2.0*TMath::Cos(para_theta) + fiberCorr_L + senDet_sizeZ/2.0);
  Tr = G4Transform3D(Ra, Ta);
  //G4VPhysicalVolume *sensitivePhysical_zp2 = new G4PVPlacement(Tr,             //Transformation
  //							      sensitive.logical, //its logical volume
  //							      "Sensitive",       //its name
  //							      world.logical,     //its mother  volume
  //							      false,	         //no boolean operation
  //							      0);	         //copy number

  
  //
  // Define quartz detector volume (1)
  //
  sec1.solid = new G4Box("Sector", sec1_sizeX/2.0, sec1_sizeY/2.0, sec1_sizeZ/2.0);
  sec1.logical = new G4LogicalVolume(sec1.solid, sec1.material, "Sector");
  Ta.setX(0.0);
  Ta.setY(0.0);
  Ta.setZ(0.0);
  Tr = G4Transform3D(Ra, Ta);
  //sec1.physical = new G4PVPlacement(Tr,            //Transformation
  //				    sec1.logical,  //its logical volume				 
  //				    "Sector",	   //its name
  //				    world.logical, //its mother  volume
  //				    false,	   //no boolean operation
  //				    0);		   //copy number
   

  //Old Paralelepiped based CpFM detector
  // Paralelepiped
  //sec1.solid = new G4Para("Sector", para_dx, para_dy, para_dz, para_alpha, para_theta, para_phi);
  /*
  sec1.solid = new G4Box("Sector", sec1_sizeX/2.0,sec1_sizeY/2.0,sec1_sizeZ/2.0);
  sec1.logical = new G4LogicalVolume(sec1.solid,sec1.material,"Sector");
  Ta.setX(0.0);
  Ta.setY(0.0);
  Ta.setZ(0.0);
  Tr = G4Transform3D(Ra, Ta);
  sec1.physical = new G4PVPlacement(Tr,            //Transformation
  				    sec1.logical,  //its logical volume				 
  				    "Sector",	     //its name
  				    world.logical, //its mother  volume
  				    false,	     //no boolean operation
				    0);      	     //copy number
  */
  //Sec4
  G4VSolid *sec4Solid = new G4Box("Sector", sec4_sizeX/2.0, sec4_sizeY/2.0, sec4_sizeZ/2.0);
  G4LogicalVolume *sec4logical = new G4LogicalVolume(sec4Solid,sec1.material,"Sector");
  G4double sec4_xinit = -(sec1_sizeZ/2.0 - sec1_sizeX)*TMath::Sin(para_theta);
  G4double sec4_yinit = 0.0;
  G4double sec4_zinit = -(sec1_sizeZ/2.0 - sec1_sizeX)*TMath::Cos(para_theta);

  sec4_xinit = sec4_xinit - sec4_sizeZ/2.0*TMath::Sin(UA9Const::angleDet);
  sec4_yinit = sec4_yinit;
  sec4_zinit = sec4_zinit - sec4_sizeZ/2.0*TMath::Cos(UA9Const::angleDet);

  sec4_xinit = sec4_xinit + (sec4_sizeX - sec1_sizeX)/2.0*TMath::Cos(UA9Const::angleDet);
  sec4_yinit = sec4_yinit;
  sec4_zinit = sec4_zinit - (sec4_sizeX - sec1_sizeX)/2.0*TMath::Sin(UA9Const::angleDet);

  Ta.setX(sec4_xinit);
  Ta.setY(sec4_yinit);
  Ta.setZ(sec4_zinit);
  Ra.rotateY(UA9Const::angleDet);
  Tr = G4Transform3D(Ra, Ta);
  if(sec4volIsOn == true){
    G4VPhysicalVolume *sec4Physical = new G4PVPlacement(Tr,            //Transformation
							sec4logical,   //its logical volume				 
							"Sector",      //its name
							world.logical, //its mother  volume
							false,	     //no boolean operation
							0);	     //copy number
  }
  Ra.rotateY(-UA9Const::angleDet);
  //1
  Ta.setX(sec4_xinit+xCopyShift);
  Ta.setY(sec4_yinit);
  Ta.setZ(sec4_zinit);
  Ra.rotateY(UA9Const::angleDet);
  Tr = G4Transform3D(Ra, Ta);
  //G4VPhysicalVolume *sec4Physical1 = new G4PVPlacement(Tr,            //Transformation
  //						       sec4logical,   //its logical volume				 
  //						       "Sector",      //its name
  //						       world.logical, //its mother  volume
  //						       false,	     //no boolean operation
  //						       0);	     //copy number
  Ra.rotateY(-UA9Const::angleDet);
  //2
  Ta.setX(sec4_xinit+2*xCopyShift);
  Ta.setY(sec4_yinit);
  Ta.setZ(sec4_zinit);
  Ra.rotateY(UA9Const::angleDet);
  Tr = G4Transform3D(Ra, Ta);
  //G4VPhysicalVolume *sec4Physical2 = new G4PVPlacement(Tr,            //Transformation
  //						       sec4logical,   //its logical volume				 
  //						       "Sector",      //its name
  //						       world.logical, //its mother  volume
  //						       false,	     //no boolean operation
  //						       0);	     //copy number
  Ra.rotateY(-UA9Const::angleDet);

  //absorber in the end of the nouse
  G4VSolid *abs2Solid = new G4Box("abs2",  abs2_sizeX/2.0,  abs2_sizeY/2.0,  abs2_sizeZ/2.0);
  G4LogicalVolume *abs2Logical = new G4LogicalVolume(abs2Solid, abs1.material, "abs2");
  //in case of the inclination angle
  G4double abs2_xinit = sec4_xinit;
  G4double abs2_yinit = sec4_yinit;
  G4double abs2_zinit = sec4_zinit;

  abs2_xinit = sec4_xinit + (sec4_sizeX + abs2_sizeX)/2.0*TMath::Cos(UA9Const::angleDet);
  abs2_yinit = sec4_yinit;
  abs2_zinit = sec4_zinit - (sec4_sizeX + abs2_sizeX)/2.0*TMath::Sin(UA9Const::angleDet);

  Ta.setX(abs2_xinit);
  Ta.setY(abs2_yinit);
  Ta.setZ(abs2_zinit);
  Ra.rotateY(UA9Const::angleDet);
  Tr = G4Transform3D(Ra, Ta);
  if(abs2volIsOn == true){
    G4VPhysicalVolume *abs2Physical = new G4PVPlacement(Tr,            //Transformation
							abs2Logical,   //its logical volume				 
							"abs2",	       //its name
							world.logical, //its mother  volume
							false,	       //no boolean operation
							0);            //copy number
  }
  Ra.rotateY(-UA9Const::angleDet);
  //1
  Ta.setX(abs2_xinit+xCopyShift);
  Ta.setY(abs2_yinit);
  Ta.setZ(abs2_zinit);
  Ra.rotateY(UA9Const::angleDet);
  Tr = G4Transform3D(Ra, Ta);
  //G4VPhysicalVolume *abs2Physical1 = new G4PVPlacement(Tr,            //Transformation
  //						      abs2Logical,  //its logical volume				 
  //						      "abs2",	  //its name
  //						      world.logical, //its mother  volume
  //						      false,	  //no boolean operation
  //						      0);		  //copy number
  Ra.rotateY(-UA9Const::angleDet);
  //2
  Ta.setX(abs2_xinit+2*xCopyShift);
  Ta.setY(abs2_yinit);
  Ta.setZ(abs2_zinit);
  Ra.rotateY(UA9Const::angleDet);
  Tr = G4Transform3D(Ra, Ta);
  //G4VPhysicalVolume *abs2Physical2 = new G4PVPlacement(Tr,            //Transformation
  //						      abs2Logical,  //its logical volume				 
  //						      "abs2",	  //its name
  //						      world.logical, //its mother  volume
  //						      false,	  //no boolean operation
  //						      0);		  //copy number
  Ra.rotateY(-UA9Const::angleDet);

  //Absorber which simulates flange 
  G4VSolid *abs3Solid = new G4Box("abs3",  abs3_sizeX/2.0,  abs3_sizeY/2.0,  abs3_sizeZ/2.0);
  G4LogicalVolume *abs3Logical = new G4LogicalVolume(abs3Solid, abs1.material, "abs3");
  G4double abs3i_xinit = sec1_sizeZ/2.0*TMath::Sin(para_theta) - sec1_sizeX/TMath::Sin(TMath::Pi()/2.0 - para_theta)/2.0;
  G4double abs3i_yinit = 0.0;
  G4double abs3i_zinit = sec1_sizeZ/2.0*TMath::Cos(para_theta);
  G4double abs3_L = 2.0*cm;//distance between flange and end of the quartz finger.
  G4double abs3_xinit = abs3i_xinit - abs3_sizeX*TMath::Cos(para_theta)/2.0 - (abs3_L + abs3_sizeZ/2.0)*TMath::Sin(para_theta);
  G4double abs3_yinit = abs3i_yinit;
  G4double abs3_zinit = abs3i_zinit + abs3_sizeX*TMath::Sin(para_theta)/2.0 - (abs3_L + abs3_sizeZ/2.0)*TMath::Cos(para_theta);
  Ta.setX(abs3_xinit);
  Ta.setY(abs3_yinit);
  Ta.setZ(abs3_zinit);
  Ra.rotateY(UA9Const::angleDet);
  Tr = G4Transform3D(Ra, Ta);
  //Left
  //G4VPhysicalVolume *abs3Physical = new G4PVPlacement(Tr,            //Transformation
  //						      abs3Logical,   //its logical volume				 
  //						      "abs3",	     //its name
  //						      world.logical, //its mother  volume
  //						      false,	     //no boolean operation
  //						      0);	     //copy number
  Ra.rotateY(-UA9Const::angleDet);
  //Right
  G4double abs3i_xinit_R = sec1_sizeZ/2.0*TMath::Sin(para_theta) + sec1_sizeX/TMath::Sin(TMath::Pi()/2.0 - para_theta)/2.0;
  G4double abs3i_yinit_R = 0.0;
  G4double abs3i_zinit_R = sec1_sizeZ/2.0*TMath::Cos(para_theta);
  G4double abs3_L_R = 2.47*cm;//distance between flange and end of the quartz finger.
  G4double abs3_xinit_R = abs3i_xinit_R + abs3_sizeX*TMath::Cos(para_theta)/2.0 - (abs3_L_R + abs3_sizeZ/2.0)*TMath::Sin(para_theta);
  G4double abs3_yinit_R = abs3i_yinit_R;
  G4double abs3_zinit_R = abs3i_zinit_R - abs3_sizeX*TMath::Sin(para_theta)/2.0 - (abs3_L_R + abs3_sizeZ/2.0)*TMath::Cos(para_theta);
  Ta.setX(abs3_xinit_R);
  Ta.setY(abs3_yinit_R);
  Ta.setZ(abs3_zinit_R);
  Ra.rotateY(UA9Const::angleDet);
  Tr = G4Transform3D(Ra, Ta);
  //Right
  //G4VPhysicalVolume *abs3Physical_R = new G4PVPlacement(Tr,            //Transformation
  //						      abs3Logical,   //its logical volume				 
  //						      "abs3",	     //its name
  //						      world.logical, //its mother  volume
  //						      false,	     //no boolean operation
  //						      0);	     //copy number
  Ra.rotateY(-UA9Const::angleDet);

  //Absorber which simulates flange 
  G4VSolid *abs4Solid = new G4Box("abs4",  abs4_sizeX/2.0,  abs4_sizeY/2.0,  abs4_sizeZ/2.0);
  G4LogicalVolume *abs4Logical = new G4LogicalVolume(abs4Solid, abs1.material, "abs4");
  G4double abs4i_xinit = sec1_sizeZ/2.0*TMath::Sin(para_theta) - sec1_sizeX/TMath::Sin(TMath::Pi()/2.0 - para_theta)/2.0;
  G4double abs4i_yinit = sec1_sizeY/2.0 + abs4_sizeX/2.0;
  G4double abs4i_zinit = sec1_sizeZ/2.0*TMath::Cos(para_theta);
  G4double abs4_L = 2.0*cm;//distance between flange and end of the quartz finger.
  G4double abs4_xinit = abs4i_xinit - abs3_sizeX*TMath::Cos(para_theta)/2.0 - (abs3_L + abs3_sizeZ/2.0)*TMath::Sin(para_theta) 
                      + (abs4_sizeY + abs3_sizeX)*TMath::Cos(para_theta)/2.0;
  G4double abs4_yinit = abs4i_yinit;
  G4double abs4_zinit = abs4i_zinit + abs3_sizeX*TMath::Sin(para_theta)/2.0 - (abs3_L + abs3_sizeZ/2.0)*TMath::Cos(para_theta)
                      - (abs4_sizeY + abs3_sizeX)*TMath::Sin(para_theta)/2.0;
  Ta.setX(abs4_xinit);
  Ta.setY(abs4_yinit);
  Ta.setZ(abs4_zinit);
  Ra.rotateZ(TMath::Pi()/2.0);
  Ra.rotateY(UA9Const::angleDet);
  Tr = G4Transform3D(Ra, Ta);
  //TOP
  //G4VPhysicalVolume *abs4Physical = new G4PVPlacement(Tr,            //Transformation
  //						      abs4Logical,   //its logical volume				 
  //						      "abs4",	     //its name
  //						      world.logical, //its mother  volume
  //						      false,	     //no boolean operation
  //						      0);	     //copy number
  Ra.rotateY(-UA9Const::angleDet);
  Ra.rotateZ(-TMath::Pi()/2.0);
  //BOTTOM
  G4double abs4i_xinit_B = sec1_sizeZ/2.0*TMath::Sin(para_theta) + sec1_sizeX/TMath::Sin(TMath::Pi()/2.0 - para_theta)/2.0;
  G4double abs4i_yinit_B = -sec1_sizeY/2.0 - abs4_sizeX/2.0;
  G4double abs4i_zinit_B = sec1_sizeZ/2.0*TMath::Cos(para_theta);
  G4double abs4_L_R = 2.47*cm;//distance between flange and end of the quartz finger.
  G4double abs4_xinit_B = abs4i_xinit_B + abs3_sizeX*TMath::Cos(para_theta)/2.0 - (abs3_L_R + abs3_sizeZ/2.0)*TMath::Sin(para_theta)
                        - (abs4_sizeY + abs3_sizeX)*TMath::Cos(para_theta)/2.0;
  G4double abs4_yinit_B = abs4i_yinit_B;
  G4double abs4_zinit_B = abs4i_zinit_B - abs3_sizeX*TMath::Sin(para_theta)/2.0 - (abs3_L_R + abs3_sizeZ/2.0)*TMath::Cos(para_theta)
                        + (abs4_sizeY + abs3_sizeX)*TMath::Sin(para_theta)/2.0;
  Ta.setX(abs4_xinit_B);
  Ta.setY(abs4_yinit_B);
  Ta.setZ(abs4_zinit_B);
  Ra.rotateZ(TMath::Pi()/2.0);
  Ra.rotateY(UA9Const::angleDet);
  Tr = G4Transform3D(Ra, Ta);
  //BOTTOM
  //G4VPhysicalVolume *abs4Physical_B = new G4PVPlacement(Tr,            //Transformation
  //						      abs4Logical,   //its logical volume				 
  //						      "abs4",	     //its name
  //						      world.logical, //its mother  volume
  //						      false,	     //no boolean operation
  //						      0);	     //copy number
  Ra.rotateY(-UA9Const::angleDet);
  Ra.rotateZ(-TMath::Pi()/2.0);


  //Sec5 for subtraction
  G4VSolid *sec5Solid = new G4Box("Sector", sec5_sizeX/2.0, sec5_sizeY/2.0, sec5_sizeZ/2.0);
  G4LogicalVolume *sec5logical = new G4LogicalVolume(sec5Solid,sec1.material,"Sector");
  //Ta.setX(-(sec1_sizeX - sec5_sizeX)/2.0);
  //Ta.setY(0.0);
  //Ta.setZ(-sec1_sizeZ/2.0 - sec5_sizeZ/2.0);
  /////////////////////////
  G4cout<<"MyCalculation"<<G4endl;
  G4RotationMatrix RaMy;
  G4ThreeVector TaMy;
  G4Transform3D TrMy;
  G4double alphaRot = twopi/4.0 - UA9Const::angleDet;
  //G4double alphaRot = 10*deg;
  TaMy.set(sec5_sizeX/2.0,0.0,sec5_sizeZ/2.0);
  //G4double zOld = TaMy.getZ() - sec1_sizeZ/2.0;
  //G4double xOld = TaMy.getX();
  G4cout<<" x "<<TaMy.getX()<<G4endl
	<<" y "<<TaMy.getY()<<G4endl
	<<" z "<<TaMy.getZ()<<G4endl;
  TaMy.rotateY(alphaRot);
  G4cout<<" x "<<TaMy.getX()<<G4endl
	<<" y "<<TaMy.getY()<<G4endl
	<<" z "<<TaMy.getZ()<<G4endl;
  G4double zNew = TaMy.getZ();
  G4double xNew = TaMy.getX();
  /////////////////////////
  Ta.setX(sec1_sizeX/2.0 - xNew);
  Ta.setY(0.0);
  Ta.setZ(-sec1_sizeZ/2.0 - zNew);
  //Ra.rotateY(twopi/4.0 - UA9Const::angleDet);
  Ra.rotateY(alphaRot);
  //Ra.rotateY(45*deg);
  Tr = G4Transform3D(Ra, Ta);
  /*
  G4VPhysicalVolume *sec5Physical = new G4PVPlacement(Tr,           //Transformation
  						      sec5logical,  //its logical volume				 
  						      "Sector",     //its name
  						      world.logical,//its mother  volume
  						      false,	    //no boolean operation
  						      0);	    //copy number
  */
  G4SubtractionSolid* subtraction =
    new G4SubtractionSolid("Sector", sec1.solid, sec5Solid,Tr);
  //G4UnionSolid* subtraction =
  //new G4UnionSolid("Sector", sec1.solid, sec5Solid,Tr);
  Ra.rotateY(-alphaRot);
  
  //Ra.rotateY(-UA9Const::angleDet);
  //Ra.rotateY(-45*deg);  

  //
  //Main sector of CpFM
  //
  G4LogicalVolume *resultlogical = new G4LogicalVolume(subtraction,sec1.material,"Sector");
  //0
  Ta.setX(0);
  Ta.setY(0);
  Ta.setZ(0);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *resultPhysical = new G4PVPlacement(Tr,            //Transformation
  							resultlogical, //its logical volume				 
  							"Sector",      //its name
  							world.logical, //its mother  volume
  							true,	       //no boolean operation
   							0);	       //copy number
  //1
  Ta.setX(xCopyShift);
  Ta.setY(0);
  Ta.setZ(zCopyShift);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *resultPhysical1 = new G4PVPlacement(Tr,            //Transformation
  							 resultlogical, //its logical volume				 
  							 "Sector",      //its name
  							 world.logical, //its mother  volume
  							 true,	        //no boolean operation
  							 0);	        //copy number
  //2
  Ta.setX(2*xCopyShift);
  Ta.setY(0);
  Ta.setZ(0);
  Tr = G4Transform3D(Ra, Ta);
  //G4VPhysicalVolume *resultPhysical2 = new G4PVPlacement(Tr,            //Transformation
  //							 resultlogical, //its logical volume				 
  //							 "Sector",      //its name
  //							 world.logical, //its mother  volume
  //							 true,	       //no boolean operation
  //							 0);	       //copy number

  abs1.solid = new G4Box("abs1",  abs1_sizeX,  abs1_sizeY,  abs1_sizeZ/2.0);
  abs1.logical = new G4LogicalVolume(abs1.solid, abs1.material, "abs1");
  //in case of the inclination angle
  Ta.setX(-sec1_sizeZ/2.0*TMath::Sin(para_theta));
  //in case of the 90 deg incl angle
  //Ta.setX(0.0);
  Ta.setY(0.0);
  //in case of the inclination angle
  Ta.setZ(-sec1_sizeZ/2.0*TMath::Cos(para_theta) - senDet_sizeZ/2.0);
  //in case of the 90 deg incl angle
  //Ta.setZ(-sec1_sizeZ/2.0 - senDet_sizeZ/2.0);
  Tr = G4Transform3D(Ra, Ta);
  //abs1.physical = new G4PVPlacement(Tr,            //Transformation
  //				    abs1.logical,  //its logical volume				 
  //				    "abs1",	   //its name
  //				    world.logical, //its mother  volume
  //				    false,	   //no boolean operation
  //				    0);		   //copy number

  /*
  //TRD 1
  G4double trd_rotAngle = (sec2_angleThetaMax - sec2_angleThetaMin)/sec2_n;
  G4double trd_dz  = sec2_d*TMath::Cos(trd_rotAngle/2.0);
  G4double trd_dx1 = 2*sec2_R*TMath::Sin(trd_rotAngle/2.0);
  G4double trd_dy1 = sec1_sizeY;  
  G4double trd_dx2 = 2*(sec2_R - sec2_d)*TMath::Sin(trd_rotAngle/2.0);
  G4double trd_mediana = (trd_dx2 + trd_dx1)/2.0; 
  G4double trd_dy2 = sec1_sizeY;
  G4double sec2_rr = (sec2_R - sec2_d)*TMath::Cos(trd_rotAngle/2.0) + trd_dz/2.0;
  //G4cout<<"trd_rotAngle "<<trd_rotAngle<<G4endl
  //	<<"trd_dz       "<<trd_dz<<G4endl
  //	<<"trd_dx1      "<<trd_dx1<<G4endl
  //	<<"trd_dy1      "<<trd_dy1<<G4endl
  //	<<"trd_dx2      "<<trd_dx2<<G4endl
  //	<<"trd_dy2      "<<trd_dy2<<G4endl;
  G4VSolid *trdSolid = new G4Trd("trdSolid", trd_dx1/2.0, trd_dx2/2.0, trd_dy1/2.0, trd_dy2/2.0, trd_dz/2.0);
  G4LogicalVolume *trdlogical = new G4LogicalVolume(trdSolid,sec1.material,"trdSolid");
  G4double trd_Angle;
  G4double trd_y0 = 0.0;
  G4double trd_z0 = -sec2_rr;
  G4double trd_x0 = 0.0;
  G4double trd_x0new;
  G4double trd_y0new = trd_y0;
  G4double trd_z0new;
  for(int i = 0;i<sec2_n;i++){
    //trd_Angle = trd_rotAngle/2.0 + trd_rotAngle*i;
    trd_Angle = trd_rotAngle/2.0 + trd_rotAngle*i;
    //Rotation
    trd_z0new = trd_z0*TMath::Cos(trd_Angle) - trd_x0*TMath::Sin(trd_Angle);
    trd_x0new = trd_z0*TMath::Sin(trd_Angle) + trd_x0*TMath::Cos(trd_Angle);
    //translation
    trd_z0new = trd_z0new + sec2_rr*TMath::Cos(trd_rotAngle/2.0) + trd_mediana/2.0*TMath::Sin(trd_rotAngle/2.0);
    Ta.setX(trd_x0new);
    Ta.setY(trd_y0new);
    Ta.setZ(trd_z0new);
    G4cout<<" trd_x0new = "<<trd_x0new<<G4endl
	  <<" trd_y0new = "<<trd_y0new<<G4endl
	  <<" trd_z0new = "<<trd_z0new<<G4endl;
    //Ra.rotateZ(90.0*deg);
    //Ra.rotateZ(-4.0*deg);
    Ra.rotateY(trd_Angle);
    Tr = G4Transform3D(Ra, Ta);
    G4VPhysicalVolume *trdPhysical = new G4PVPlacement(Tr,            //Transformation
						       trdlogical,    //its logical volume				 
						       "trdSolid",    //its name
						       world.logical, //its mother  volume
						       false,	      //no boolean operation
						       0);	      //copy number
    Ra.rotateY(-trd_Angle);
    //Ra.rotateZ(4.0*deg);
    //Ra.rotateZ(-90.0*deg);
  }
  */  

  /////////FIBERS///////////

  fiberCorr.solid = new G4Tubs("Sector", fiberCorr_Rmin, fiberCorr_Rmax, fiberCorr_L/2.0, 0.0,  360.0*deg);
  fiberCorr.logical = new G4LogicalVolume(fiberCorr.solid, fiberCorr.material, "Sector");
  fiberClad.solid = new G4Tubs("SectorClad", fiberClad_Rmin, fiberClad_Rmax, fiberClad_L/2.0, 0.0,  360.0*deg);
  fiberClad.logical = new G4LogicalVolume(fiberClad.solid, fiberClad.material, "SectorClad");
  fiberCoat.solid = new G4Tubs("SectorCoat", fiberCoat_Rmin, fiberCoat_Rmax, fiberCoat_L/2.0, 0.0,  360.0*deg);
  fiberCoat.logical = new G4LogicalVolume(fiberCoat.solid, fiberCoat.material, "SectorCoat");

  
  int ii = 0;//along z dirrection
  int jj = 0;//along y dirrection 

  G4double fiber_X0 = sec1_sizeZ/2.0*TMath::Sin(para_theta) - senDet_sizeX + fiberCoat_Rmax;
  G4double fiber_Y0 = -sec1_sizeY/2.0 + fiberCoat_Rmax;
  G4double fiber_Z0 = sec1_sizeZ/2.0*TMath::Cos(para_theta) + fiberCorr_L/2.0;
  G4double fiber_dd = 0.02*mm;
  //0
  for(jj = 0;jj<16;jj++){
    //for(ii = 0;ii<19;ii++){
    for(ii = 0;ii<10;ii++){
      if(jj%2 == 0){
	Ta.setX(fiber_X0 + ii*2*(fiberCoat_Rmax + fiber_dd));
      }
      else{
	Ta.setX(fiber_X0 + ii*2*(fiberCoat_Rmax + fiber_dd) - 2*(fiberCoat_Rmax + fiber_dd)*TMath::Sin(30.0*deg));
      }
      Ta.setY(fiber_Y0 + jj*2*(fiberCoat_Rmax + fiber_dd)*TMath::Cos(30.0*deg));
      Ta.setZ(fiber_Z0);
      Tr = G4Transform3D(Ra, Ta);
      //fiberCorr.physical = new G4PVPlacement( Tr, fiberCorr.logical, "Sector", world.logical, false, 0);
      //fiberClad.physical = new G4PVPlacement( Tr, fiberClad.logical, "SectorClad", world.logical, false, 0);
      //fiberCoat.physical = new G4PVPlacement( Tr, fiberCoat.logical, "SectorCoat", world.logical, false, 0);
    }
  }
  //1
  for(jj = 0;jj<16;jj++){
    for(ii = 0;ii<6;ii++){
      if(jj%2 == 0){
	Ta.setX(fiber_X0 + xCopyShift + ii*2*(fiberCoat_Rmax + fiber_dd));
      }
      else{
	Ta.setX(fiber_X0 + xCopyShift + ii*2*(fiberCoat_Rmax + fiber_dd) - 2*(fiberCoat_Rmax + fiber_dd)*TMath::Sin(30.0*deg));
      }
      Ta.setY(fiber_Y0 + jj*2*(fiberCoat_Rmax + fiber_dd)*TMath::Cos(30.0*deg));
      Ta.setZ(fiber_Z0);
      Tr = G4Transform3D(Ra, Ta);
      //fiberCorr.physical = new G4PVPlacement( Tr, fiberCorr.logical, "Sector", world.logical, false, 0);
      //fiberClad.physical = new G4PVPlacement( Tr, fiberClad.logical, "SectorClad", world.logical, false, 0);
      //fiberCoat.physical = new G4PVPlacement( Tr, fiberCoat.logical, "SectorCoat", world.logical, false, 0);
    }
  }
  //2
  for(jj = 0;jj<16;jj++){
    for(ii = 0;ii<6;ii++){
      if(jj%2 == 0){
	Ta.setX(fiber_X0 + 2*xCopyShift + ii*2*(fiberCoat_Rmax + fiber_dd));
      }
      else{
	Ta.setX(fiber_X0 + 2*xCopyShift + ii*2*(fiberCoat_Rmax + fiber_dd) - 2*(fiberCoat_Rmax + fiber_dd)*TMath::Sin(30.0*deg));
      }
      Ta.setY(fiber_Y0 + jj*2*(fiberCoat_Rmax + fiber_dd)*TMath::Cos(30.0*deg));
      Ta.setZ(fiber_Z0);
      Tr = G4Transform3D(Ra, Ta);
      //fiberCorr.physical = new G4PVPlacement( Tr, fiberCorr.logical, "Sector", world.logical, false, 0);
      //fiberClad.physical = new G4PVPlacement( Tr, fiberClad.logical, "SectorClad", world.logical, false, 0);
      //fiberCoat.physical = new G4PVPlacement( Tr, fiberCoat.logical, "SectorCoat", world.logical, false, 0);
    }
  }


  //
  // Set Visualization Attributes
  //
  G4Color blue        = G4Color(0., 0., 1.);
  G4Color green       = G4Color(0., 1., 0.);
  G4Color red         = G4Color(1., 0., 0.);
  G4Color white       = G4Color(1., 1., 1.);
  G4Color cyan        = G4Color(0., 1., 1.);
  G4Color DircColor   = G4Color(0.0, 0.0, 1.0, 0.2);
  G4Color SensColor   = G4Color(0.0, 1.0, 1.0, 0.1);

  worldVisAtt->SetColor(white);
  worldVisAtt->SetVisibility(true);
  quartzVisAtt->SetColor(DircColor);
  quartzVisAtt->SetVisibility(true);
  sensitiveVisAtt->SetColor(SensColor);
  sensitiveVisAtt->SetVisibility(true);
  absVisAtt->SetColor(red);
  absVisAtt->SetVisibility(true);

  //trdlogical->SetVisAttributes(quartzVisAtt);
  sec1.logical->SetVisAttributes(quartzVisAtt);
  sec4logical->SetVisAttributes(quartzVisAtt);
  resultlogical->SetVisAttributes(quartzVisAtt);
  sec5logical->SetVisAttributes(quartzVisAtt);

  sensitive.logical->SetVisAttributes(sensitiveVisAtt);
  abs1.logical->SetVisAttributes(absVisAtt);
  abs2Logical->SetVisAttributes(absVisAtt);
  abs3Logical->SetVisAttributes(absVisAtt);
  abs4Logical->SetVisAttributes(absVisAtt);

  fiberCorr.logical->SetVisAttributes(quartzVisAtt);
  fiberClad.logical->SetVisAttributes(quartzVisAtt);
  fiberCoat.logical->SetVisAttributes(absVisAtt);

  //world.logical->SetVisAttributes(worldVisAtt);
  
  //
  // Define Optical Borders
  //

  // Surface for killing photons at borders
  const G4int num1 = 2;
  G4double Ephoton[num1] = {1.5*eV, 5.8*eV};

  G4OpticalSurface* OpVolumeKillSurface =
    new G4OpticalSurface("VolumeKillSurface");
  OpVolumeKillSurface->SetType(dielectric_metal);
  OpVolumeKillSurface->SetFinish(polished);
  OpVolumeKillSurface->SetModel(glisur);


  G4double ReflectivityKill[num1] = {0., 0.};
  G4double EfficiencyKill[num1] = {1., 1.};
  G4MaterialPropertiesTable* VolumeKill = new G4MaterialPropertiesTable();
  VolumeKill->AddProperty("REFLECTIVITY", Ephoton, ReflectivityKill, num1);
  VolumeKill->AddProperty("EFFICIENCY",   Ephoton, EfficiencyKill,   num1);
  OpVolumeKillSurface->SetMaterialPropertiesTable(VolumeKill);
  new G4LogicalSkinSurface("SensitiveSurface", 
			   sensitive.logical, OpVolumeKillSurface);
  
  
  // 
  // Sensitive detector definition
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  SensitiveDetector* aSD = new SensitiveDetector("fTOF");
  SDman->AddNewDetector(aSD);
  //sensitive.logical->SetSensitiveDetector(aSD);
  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  G4double maxStep   = 0.1*mm;
  G4double maxLength = 5.0*m;
  G4double maxTime   = 20.0*ns; 
  G4double minEkin   = 1.0/100*MeV;
  G4double mionRang  = 0.01*mm;
  stepLimit = new G4UserLimits(maxStep,maxLength,maxTime,minEkin,mionRang);
  sec1.logical->SetUserLimits(stepLimit);

  //G4GDMLParser parser;
  //parser.Write("CpFM.gdml", world.physical);
  //parser.Write("CpFM_t.gdml", world.physical);

  return world.physical;
}
