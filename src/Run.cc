#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4ProcessTable.hh"
#include "G4Radioactivation.hh"
#include "G4TwoVector.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det), fParticle(0), fEkin(0.)
{
  fEdep = 0.;
}

Run::~Run()
{}

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
}

void Run::ParticleCount(G4String name, G4double Ekin)
{
  std::map<G4String, ParticleData>::iterator it = fParticleDataMap.find(name);
  if ( it == fParticleDataMap.end()) {
    fParticleDataMap[name] = ParticleData(1, Ekin, Ekin, Ekin);
  }
   else {
     ParticleData& data = it->second;
     data.fCount++;
     data.fEmean += Ekin;
     //update min max
     G4double emin = data.fEmin;
     if (Ekin < emin) data.fEmin = Ekin;
     G4double emax = data.fEmax;
     if (Ekin > emax) data.fEmax = Ekin; 
   }   
 }
 
void Run::AddEdep(G4double edep)
{ 
  fEdep += edep;
}

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  //primary particle info
  //
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;
    
  //map: created particles  
  std::map<G4String,ParticleData>::const_iterator itc;
  for (itc = localRun->fParticleDataMap.begin(); 
       itc != localRun->fParticleDataMap.end(); ++itc) {
    
    G4String name = itc->first;
    const ParticleData& localData = itc->second;   
    if ( fParticleDataMap.find(name) == fParticleDataMap.end()) {
      fParticleDataMap[name]
       = ParticleData(localData.fCount, 
                      localData.fEmean, 
                      localData.fEmin, 
                      localData.fEmax);
    }
    else {
      ParticleData& data = fParticleDataMap[name];   
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin) data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax) data.fEmax = emax; 
    }   
  }

  G4Run::Merge(run); 
} 

void Run::EndOfRun() 
{
  G4int prec = 5, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);
  
  // run condition
  //   
  G4String Particle = fParticle->GetParticleName();    
  G4cout << "\n The run is " << numberOfEvent << " "<< Particle << " of "
         << G4BestUnit(fEkin,"Energy") << " through : ";

  G4cout << G4endl;
    
  // particles count
  //
  G4cout << "\n List of generated particles:" << G4endl;
     
  std::map<G4String,ParticleData>::iterator itc;               
  for (itc = fParticleDataMap.begin(); itc != fParticleDataMap.end(); itc++) {
     G4String name = itc->first;
     ParticleData data = itc->second;
     G4int count = data.fCount;
     G4double eMean = data.fEmean/count;
     G4double eMin = data.fEmin;
     G4double eMax = data.fEmax;    
         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ")" << G4endl;         
 }
  G4cout << G4endl;
 
  // activities in VR mode
  //
  WriteActivity(numberOfEvent);
  fParticleDataMap.clear();    
                          
  //restore default format         
  G4cout.precision(dfprec);   
}

void Run::WriteActivity(G4int nevent)
{
 G4ProcessTable *pTable = G4ProcessTable::GetProcessTable();
 G4Radioactivation* rDecay = (G4Radioactivation *)
         pTable->FindProcess("Radioactivation", "GenericIon");
   
 // output the induced radioactivities (in VR mode only)
 //
 if ((rDecay == 0) || (rDecay->IsAnalogueMonteCarlo())) return;
 
 G4String fileName = G4AnalysisManager::Instance()->GetFileName() + ".activity";
 std::ofstream outfile (fileName, std::ios::out );
 
 std::vector<G4RadioactivityTable*> theTables =
                              rDecay->GetTheRadioactivityTables();

 for (size_t i = 0 ; i < theTables.size(); i++) {
    G4double rate, error;
    outfile << "Radioactivities in decay window no. " << i << G4endl;
    outfile << "Z \tA \tE \tActivity (decays/window) \tError (decays/window) "
            << G4endl;

    map<G4ThreeVector,G4TwoVector> *aMap = theTables[i]->GetTheMap();
    map<G4ThreeVector,G4TwoVector>::iterator iter;
    for (iter=aMap->begin(); iter != aMap->end(); iter++) {
       rate = iter->second.x()/nevent;
       error = std::sqrt(iter->second.y())/nevent;
       if (rate < 0.) rate = 0.;                // statically it can be < 0.
       outfile << iter->first.x() <<"\t"<< iter->first.y() <<"\t"
               << iter->first.z() << "\t" << rate <<"\t" << error << G4endl;
    }
    outfile << G4endl;
 }
 outfile.close();
}
