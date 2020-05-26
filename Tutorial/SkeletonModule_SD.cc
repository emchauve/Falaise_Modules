#include <bayeux/dpp/chain_module.h>
#include <bayeux/mctools/simulated_data.h>

class SkeletonModule_SD : public dpp::chain_module
{
public:
  SkeletonModule_SD();
  ~SkeletonModule_SD();
  
  void initialize(const datatools::properties &, datatools::service_manager &, dpp::module_handle_dict_type &);
  dpp::chain_module::process_status process(datatools::things &);

  void finalize();
  
private:
  DPP_MODULE_REGISTRATION_INTERFACE(SkeletonModule_SD);

};

DPP_MODULE_REGISTRATION_IMPLEMENT(SkeletonModule_SD, "Skeleton_SD");


SkeletonModule_SD::SkeletonModule_SD() : dpp::chain_module() {}

SkeletonModule_SD::~SkeletonModule_SD() {this->finalize();}

void SkeletonModule_SD::initialize(const datatools::properties &, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  std::cout << "+++ SkeletonModule_SD::initialize()" << std::endl;

  this->_set_initialized(true);
}

dpp::chain_module::process_status SkeletonModule_SD::process (datatools::things &event)
{
  std::cout << "+++ SkeletonModule_SD::process()" << std::endl;

  const mctools::simulated_data & sd_bank = event.get<mctools::simulated_data>("SD");
  // sd_bank.tree_dump(std::cout);

  ///////////////////
  // PRIMARY EVENT //
  ///////////////////
  
  const mctools::simulated_data::primary_event_type & primary_event = sd_bank.get_primary_event();

  // primary particle
  const genbb::primary_event::particles_col_type & primary_particles = primary_event.get_particles();

  int nb_electrons = 0;
  double total_energy = 0;
  
  for (const genbb::primary_particle & a_primary_particle : primary_particles)
    {
      // a_primary_particle.get_type();
      // a_primary_particle.get_pdg_code();
      // a_primary_particle.get_kinetic_energy();

      if (a_primary_particle.get_type() == genbb::primary_particle::ELECTRON)
	nb_electrons++;
      
      total_energy += a_primary_particle.get_kinetic_energy();
    }

  std::cout << "+++ SkeletonModule_SD::process() : nb_electrons = " << nb_electrons << std::endl;
  std::cout << "+++ SkeletonModule_SD::process() : total_energy = " << total_energy << std::endl;

  // primary_event.get_genbb_weight();

  ////////////////////
  // PRIMARY VERTEX //
  ////////////////////
  
  const geomtools::vector_3d & primary_vertex = sd_bank.get_vertex();

  // primary_vertex.x()/CLHEP::mm;
  // primary_vertex.y()/CLHEP::mm;
  // primary_vertex.z()/CLHEP::mm;

  //////////////////////////////////
  // CALORIMETER+TRACKER STEP HIT //
  //////////////////////////////////
  
  const mctools::simulated_data::step_hits_dict_type & primary_step_hits = sd_bank.get_step_hits_dict();

  if (sd_bank.has_step_hits("gg"))
    {
      mctools::simulated_data::hit_handle_collection_type step_hits = sd_bank.get_step_hits("gg");

      int step_hit_index = 0;
      for (const mctools::simulated_data::hit_handle_type a_step_hit : step_hits)
	{
	  std::cout << std::endl << "+++ step hit \"gg\" #" << step_hit_index << std::endl;
	  a_step_hit.get().dump();
	}
      
      // const mctools::simulated_data::base_step_hit 
    }
  
  return PROCESS_OK;
}

void SkeletonModule_SD::finalize()
{
  std::cout << "+++ SkeletonModule_SD::finalize()" << std::endl;
  this->_set_initialized(false);
}
