#include <bayeux/dpp/chain_module.h>
#include <bayeux/mctools/simulated_data.h>

class FLMANU_Test : public dpp::chain_module
{
public:
  FLMANU_Test();
  virtual ~FLMANU_Test();

  virtual void initialize(const datatools::properties &myConfig,
                          datatools::service_manager& flServices,
                          dpp::module_handle_dict_type& moduleDict);

  virtual dpp::chain_module::process_status process(datatools::things& workItem);

  virtual void finalize();
  
private:
  DPP_MODULE_REGISTRATION_INTERFACE(FLMANU_Test);

  unsigned long long total_entries;
  unsigned long long select_entries;
};

DPP_MODULE_REGISTRATION_IMPLEMENT(FLMANU_Test, "FLMANU_Test");

FLMANU_Test::FLMANU_Test() : dpp::chain_module() {}

FLMANU_Test::~FLMANU_Test() {this->finalize();}

void FLMANU_Test::initialize(const datatools::properties  &myConfig,
			      datatools::service_manager   &/*&flServices*/,
			      dpp::module_handle_dict_type &/*moduleDict*/)
{
  total_entries = 0; 
  select_entries = 0;

  this->_set_initialized(true);
}

dpp::chain_module::process_status FLMANU_Test::process (datatools::things &event)
{
  ++total_entries;

  const mctools::simulated_data &simulated_data = event.get<mctools::simulated_data> ("SD");

  int calorimeter_hit = 0;
  
  if (simulated_data.has_step_hits("calo"))
    calorimeter_hit += simulated_data.get_number_of_step_hits("calo");
    
  if (simulated_data.has_step_hits("xcalo"))
    calorimeter_hit += simulated_data.get_number_of_step_hits("xcalo");
    
  if (simulated_data.has_step_hits("gveto"))
    calorimeter_hit += simulated_data.get_number_of_step_hits("gveto");
  
  if (calorimeter_hit > 0) PROCESS_OK;

  else return PROCESS_OK;
  
}

void FLMANU_Test::finalize()
{
  std::cout << "FLMANU_Test::finalize() : " << select_entries << "/" << total_entries << " selected" << std::endl;

  this->_set_initialized(false);
}
