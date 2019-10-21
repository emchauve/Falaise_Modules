#include <bayeux/dpp/chain_module.h>
#include <bayeux/mctools/simulated_data.h>

class zero_hit_cleaner : public dpp::chain_module
{
public:
  zero_hit_cleaner ()
  {
    // default cut value
    calorimeter_hit_threshold = 1;
    tracker_hit_threshold = 1;
  }

  virtual void initialize (const datatools::properties & properties,
                           datatools::service_manager &, dpp::module_handle_dict_type &)
  {
    if (properties.has_key("calorimeter_hit_threshold"))
      calorimeter_hit_threshold = properties.fetch_integer("calorimeter_hit_threshold");

    if (properties.has_key("tracker_hit_threshold"))
      tracker_hit_threshold = properties.fetch_integer("tracker_hit_threshold");
    
    total_events = selected_events = 0;

    std::cout << "ZeroHitCleaner : gg_hit >= " << tracker_hit_threshold
	      << " and calo_hit >= " << calorimeter_hit_threshold << std::endl;
  }

  ~zero_hit_cleaner()
  {
    if (total_events > 0)
      std::cout << "zeroHitCleaner summary : "  << selected_events << " / " << total_events << " events selected"
		<< " (" <<  (100.0*selected_events)/total_events << " %)" << std::endl;
  }

  dpp::chain_module::process_status process (datatools::things &event)
  {
    ++total_events;
    
    auto & simulated_data = event.get<mctools::simulated_data>("SD");

    // calorimeter nhit cut //
    
    int calorimeter_hit = 0;

    if (simulated_data.has_step_hits("calo"))
      calorimeter_hit += simulated_data.get_number_of_step_hits("calo");

    if (simulated_data.has_step_hits("xcalo"))
      calorimeter_hit += simulated_data.get_number_of_step_hits("xcalo");

    if (simulated_data.has_step_hits("gveto"))
      calorimeter_hit += simulated_data.get_number_of_step_hits("gveto");
    
    if (calorimeter_hit < calorimeter_hit_threshold)
      return PROCESS_STOP;

    // tracker nhit cut //
    
    int tracker_hit = 0;

    if (simulated_data.has_step_hits("gg"))
      tracker_hit += simulated_data.get_number_of_step_hits("gg");
    
    if (tracker_hit < tracker_hit_threshold)
      return PROCESS_STOP;

    // event passed nhit cut //
    
    ++selected_events;
    
    return PROCESS_OK;
  }

private:
  DPP_MODULE_REGISTRATION_INTERFACE(zero_hit_cleaner);

  int tracker_hit_threshold;
  int calorimeter_hit_threshold;
  
  int total_events;
  int selected_events;
};

DPP_MODULE_REGISTRATION_IMPLEMENT(zero_hit_cleaner, "zero_hit_cleaner")
