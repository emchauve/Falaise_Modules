#include <bayeux/dpp/chain_module.h>

#include <bayeux/mctools/simulated_data.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>

class FLMANU_GENBB : public dpp::chain_module
{
public:
  FLMANU_GENBB()
  {
    binSpectrum = 500;

    minSpectrum = 0.0;
    maxSpectrum = 5.0;

    output_filename = "flmanu-genbb-output.root";
    
    this->_set_initialized(false);
  }
  
  ~FLMANU_GENBB()
  {
    if (this->is_initialized())
      this->finalize();
  }
    
  virtual void initialize (const datatools::properties & properties,
                           datatools::service_manager &, dpp::module_handle_dict_type &)
  {
    if (properties.has_key("nbins"))
      binSpectrum = properties.fetch_integer("nbins");

    if (properties.has_key("hmin"))
      minSpectrum = properties.fetch_real("hmin");

    if (properties.has_key("hmax"))
      maxSpectrum = properties.fetch_real("hmax");

    if (properties.has_key("output_filename"))
      output_filename = properties.fetch_string("output_filename");
  
    std::cout << "FLMANU_GENBB : output file = " << output_filename << std::endl;
    
    outputFile = new TFile (output_filename.data(), "RECREATE");

    // ROOT tree

    outputTree = new TTree ("genbb", "");

    outputTree->Branch("weight", &outputEventWeight);

    outputTree->Branch("x", &outputVertexX);
    outputTree->Branch("y", &outputVertexY);
    outputTree->Branch("z", &outputVertexZ);

    outputParticleType = new std::vector<unsigned short>;
    outputTree->Branch("type",   &outputParticleType);

    outputParticleTime  = new std::vector<unsigned short>;
    outputTree->Branch("time",   &outputParticleTime);

    outputParticlePx = new std::vector<float>;
    outputParticlePy = new std::vector<float>;
    outputParticlePz = new std::vector<float>;
    outputTree->Branch("px", &outputParticlePx);
    outputTree->Branch("py", &outputParticlePy);
    outputTree->Branch("pz", &outputParticlePz);

    outputParticleEnergy = new std::vector<float>;
    outputTree->Branch("energy", &outputParticleEnergy);

    // ROOT histograms 

    weightSpectrum   = new TH1F ("weight_spectrum", "event weight", 1000, 0, 1);
      
    gammaSpectrum    = new TH1F ("gamma_energy", "#gamma spectrum;Energy (MeV);",     binSpectrum, minSpectrum, maxSpectrum);
    electronSpectrum = new TH1F ("electron_energy", "#beta^{-} spectrum;Energy (MeV);", binSpectrum, minSpectrum, maxSpectrum);
    positronSpectrum = new TH1F ("positron_energy", "#beta^{+} spectrum;Energy (MeV);", binSpectrum, minSpectrum, maxSpectrum);
    alphaSpectrum    = new TH1F ("alpha_energy", "#alpha spectrum;Energy (MeV);",     binSpectrum, minSpectrum, maxSpectrum);
    totalSpectrum    = new TH1F ("total_energy", "#gamma #spectrum;Energy (MeV);",    binSpectrum, minSpectrum, maxSpectrum);

    // xyVertex = new TH2F ("xy_vertex", "XY vertex", 2000, -4, 4, 2000, -4, 4);
    // xzVertex = new TH2F ("xz_vertex", "XZ vertex", 2000, -4, 4, 2000, -4, 4);
    // yzVertex = new TH2F ("yz_vertex", "yz vertex", 2000, -4, 4, 2000, -4, 4);

    this->_set_initialized(true);
  }
 
  void finalize ()
  {
    outputFile->cd();

    outputTree->Write();

    weightSpectrum->Write();
    
    gammaSpectrum->Write();
    electronSpectrum->Write();
    positronSpectrum->Write();
    alphaSpectrum->Write();
    totalSpectrum->Write();

    // xyVertex->Write();
    // xzVertex->Write();
    // yzVertex->Write();
    
    outputFile->Close();    
  }
    

  dpp::chain_module::process_status process(datatools::things &event)
  {     
    auto & simulatedData = event.get<mctools::simulated_data>("SD");

    auto & primaryEvent = simulatedData.get_primary_event();

    // GENBB weight
    
    outputEventWeight = primaryEvent.get_genbb_weight();
    weightSpectrum->Fill(primaryEvent.get_genbb_weight());
    
    // GENBB vertex

    const geomtools::vector_3d & primaryVertex = simulatedData.get_vertex();

    outputVertexX = primaryVertex.x()/CLHEP::mm;
    outputVertexY = primaryVertex.y()/CLHEP::mm;
    outputVertexZ = primaryVertex.z()/CLHEP::mm;

    // GENBB particles

    outputParticleType->clear();
    outputParticleTime->clear();
    outputParticlePx->clear();
    outputParticlePy->clear();
    outputParticlePz->clear();
    outputParticleEnergy->clear();

    float totalKineticEnergy = 0;
    
    auto & primaryParticles =  	primaryEvent.get_particles();

    for (const auto & primaryParticle : primaryParticles)
      {
    	float kineticEnergy = primaryParticle.get_kinetic_energy();
    	totalKineticEnergy += kineticEnergy;

        outputParticleType->push_back   ( primaryParticle.get_type() );
	outputParticleTime->push_back   ( primaryParticle.get_time() );
        outputParticlePx->push_back     ( primaryParticle.get_momentum()[0] );
        outputParticlePy->push_back     ( primaryParticle.get_momentum()[1] );
	outputParticlePz->push_back     ( primaryParticle.get_momentum()[2] );
        outputParticleEnergy->push_back ( primaryParticle.get_kinetic_energy() );

    	switch (primaryParticle.get_type())
    	  {
    	  case genbb::primary_particle::ELECTRON :
    	    electronSpectrum->Fill(kineticEnergy);
    	    break;

    	  case genbb::primary_particle::GAMMA :
    	    gammaSpectrum->Fill(kineticEnergy);
    	    break;

    	  case genbb::primary_particle::POSITRON :
    	    positronSpectrum->Fill(kineticEnergy);
    	    break;

    	  case genbb::primary_particle::ALPHA :
    	    alphaSpectrum->Fill(kineticEnergy);
    	    break;
	    
    	  default:
    	    break;
    	  }
      }

    totalSpectrum->Fill(totalKineticEnergy);
    
    outputTree->Fill();

    // // vertex

    // const geomtools::vector_3d & primaryVertex = simulatedData.get_vertex();

    // float vertexX = primaryVertex.x()/CLHEP::m;
    // float vertexY = primaryVertex.y()/CLHEP::m;
    // float vertexZ = primaryVertex.z()/CLHEP::m;
      
    // xyVertex->Fill(vertexX, vertexY);
    // xzVertex->Fill(vertexX, vertexZ);
    // yzVertex->Fill(vertexY, vertexZ);
    
    return PROCESS_OK;
  }

private:
  DPP_MODULE_REGISTRATION_INTERFACE(FLMANU_GENBB);

  bool isInitialized;
    
  std::string output_filename;

  int   binSpectrum;
  float minSpectrum;
  float maxSpectrum;
  
  TFile *outputFile;

  TTree *outputTree;

  float outputVertexX;
  float outputVertexY;
  float outputVertexZ;

  float outputEventWeight;

  std::vector<unsigned short> *outputParticleType;
  std::vector<unsigned short> *outputParticleTime;
  std::vector<float>          *outputParticlePx;
  std::vector<float>          *outputParticlePy;
  std::vector<float>          *outputParticlePz;
  std::vector<float>          *outputParticleEnergy;

  TH1F *weightSpectrum;
  
  TH1F *gammaSpectrum;
  TH1F *electronSpectrum;
  TH1F *positronSpectrum;
  TH1F *alphaSpectrum;
  TH1F *totalSpectrum;

  TH2F *xyVertex;
  TH2F *xzVertex;
  TH2F *yzVertex; 
};

DPP_MODULE_REGISTRATION_IMPLEMENT(FLMANU_GENBB, "FLMANU_GENBB")
