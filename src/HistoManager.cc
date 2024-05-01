#include "HistoManager.hh"
#include "G4UnitsTable.hh"

HistoManager::HistoManager()
  : fFileName("rdecay02")
{
  Book();
}

HistoManager::~HistoManager()
{
}

void HistoManager::Book()
{
  // Create or get analysis manager
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  
  analysis->SetDefaultFileType("root");
  analysis->SetFileName(fFileName);
  analysis->SetVerboseLevel(1);
  analysis->SetActivation(true);     //enable inactivation of histos, nTuples
    
  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  //
  ////analysis->SetHistoDirectoryName("histo");  
  ////analysis->SetFirstHistoId(1);
    
  G4int id = analysis->CreateH1("H10","Energy deposit (MeV)",
                       nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);
    
  /*id = analysis->CreateH1("H11","Energy deposit (MeV) in the detector",
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);*/

  id = analysis->CreateH1("H12","Total energy (MeV)",
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);  

  id = analysis->CreateH1("H16","Decay emission spectrum (0 - 10 MeV)",
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);  
  
  id = analysis->CreateH1("H17","Decay emission spectrum (0 - 1 MeV)",
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);

  id = analysis->CreateH1("H18","Decay emission spectrum (0 - 0.1 MeV)",
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);
  
  // nTuple
  analysis->CreateNtuple("RDecayProducts", "All Products of RDecay");
  analysis->CreateNtupleDColumn("PID");       //column 0
  analysis->CreateNtupleDColumn("Z");         //column 1
  analysis->CreateNtupleDColumn("A");         //column 2    
  analysis->CreateNtupleDColumn("Energy");    //column 3
  analysis->FinishNtuple();
  
  analysis->SetNtupleActivation(true);          
}
