#include <bayeux/dpp/chain_module.h>

//////////////////////////////
// Module class declaration //
//////////////////////////////

class SkeletonModule : public dpp::chain_module
{
public:
  SkeletonModule();
  ~SkeletonModule();
  
  void initialize(const datatools::properties &, datatools::service_manager &, dpp::module_handle_dict_type &);
  dpp::chain_module::process_status process(datatools::things &);

  void finalize();
  
private:
  // Bayeux' macro to register this class as data processing module
  DPP_MODULE_REGISTRATION_INTERFACE(SkeletonModule);

  unsigned long long nb_events_processed;
};

/////////////////////////////
// Module class definition //
/////////////////////////////

// Bayeux' macro to register this class as data processing module
DPP_MODULE_REGISTRATION_IMPLEMENT(SkeletonModule, "Skeleton");


SkeletonModule::SkeletonModule() : dpp::chain_module() {}

SkeletonModule::~SkeletonModule()
{
  if (this->is_initialized())
    this->finalize();
}

void SkeletonModule::initialize(const datatools::properties &, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  std::cout << "+++ SkeletonModule::initialize()" << std::endl;
  
  nb_events_processed = 0; 
  this->_set_initialized(true);
}

dpp::chain_module::process_status SkeletonModule::process (datatools::things &event)
{
  std::cout << "+++ SkeletonModule::process()" << std::endl;
  ++nb_events_processed;
  return PROCESS_OK;
  
}

void SkeletonModule::finalize()
{
  std::cout << "+++ SkeletonModule::finalize()" << std::endl;
  this->_set_initialized(false);
}
