#include <bayeux/dpp/chain_module.h>
#include <bayeux/genbb_help/primary_event.h>
#include <bayeux/mctools/base_step_hit.h>
#include <bayeux/mctools/simulated_data.h>

#include <falaise/snemo/datamodels/event_header.h>
#include <falaise/snemo/datamodels/calibrated_data.h>

class FLMANU_Killer : public dpp::chain_module
{
public:
  FLMANU_Killer();
  virtual ~FLMANU_Killer();
  virtual void initialize(const datatools::properties &myConfig,
                          datatools::service_manager& flServices,
                          dpp::module_handle_dict_type& moduleDict);

  virtual dpp::chain_module::process_status process(datatools::things &event);
  
  virtual void finalize ();

private:
  DPP_MODULE_REGISTRATION_INTERFACE(FLMANU_Killer);

  bool om_num_to_kill[712];
  bool cell_num_to_kill[2034];
};

DPP_MODULE_REGISTRATION_IMPLEMENT(FLMANU_Killer, "FLMANU_Killer");

FLMANU_Killer::FLMANU_Killer() : dpp::chain_module() {}

FLMANU_Killer::~FLMANU_Killer()
{
  this->finalize();
}

void FLMANU_Killer::initialize(const datatools::properties  &myConfig,
			      datatools::service_manager   &/*&flServices*/,
			      dpp::module_handle_dict_type &/*moduleDict*/)
{
  if (myConfig.has_key("kill_all_cd"))
    {
      bool kill = myConfig.fetch_boolean("kill_all_cd");
      DT_LOG_DEBUG(get_logging_priority(), "setting kill all CD to " << kill);
      memset(om_num_to_kill, kill, 712*sizeof(bool));
      memset(cell_num_to_kill, kill, 2034*sizeof(bool));
    }

  if (myConfig.has_key("keep_area0_cd"))
    {
      if (myConfig.fetch_boolean("keep_area0_cd"))
	{
	  // enable cells
	  for (int side=0; side<2; ++side)
	    for (int row=0; row<14; ++row)
	      for (int layer=0; layer<9; ++layer)
		cell_num_to_kill[1017*side+9*row+layer] = false;

	  // enable MW OMs
	  for (int side=0; side<2; ++side)
	    for (int column=0; column<3; ++column)
	      for (int row=0; row<13; ++row)
		om_num_to_kill[side*260+column*13+row] = false;

	  // enable XW OMs on mountain side
	  for (int side=0; side<2; ++side) // side = 0,1
	    for (int wall=0; wall<1; ++wall) // wall 0 only
	    for (int column=0; column<2; ++column)
	      for (int row=0; row<16; ++row)
		om_num_to_kill[520+side*64+wall*32+column*16+row] = false;
	}
    }
    
  this->_set_initialized(true);
}

dpp::chain_module::process_status FLMANU_Killer::process(datatools::things &event)
{
  snemo::datamodel::calibrated_data & CD = event.grab<snemo::datamodel::calibrated_data>("CD");
  
  // for (const auto & tracker_hit : CD.tracker_hits())
  for (auto tracker_hit = CD.tracker_hits().begin(); tracker_hit != CD.tracker_hits().end();)
    {
      unsigned short a_cell_num = 9 * 113 * (*tracker_hit)->get_side();
      a_cell_num += 9 * (*tracker_hit)->get_row();
      a_cell_num += (*tracker_hit)->get_layer();
 
      if (cell_num_to_kill[a_cell_num])	
	CD.tracker_hits().erase(tracker_hit);
      else ++tracker_hit;
    }


  // for (const auto & calo_hit : CD.calorimeter_hits())
  for (auto calo_hit = CD.calorimeter_hits().begin(); calo_hit != CD.calorimeter_hits().end();)
    {
      unsigned short an_om_num = -1;

      const geomtools::geom_id & calo_geomid = (*calo_hit)->get_geom_id();

      if (calo_geomid.get_type() == 1302) // mwall
	an_om_num = calo_geomid.get(1)*20*13 + calo_geomid.get(2)*13 + calo_geomid.get(3);
      else if (calo_geomid.get_type() == 1232) // xwall
	an_om_num = 520 + calo_geomid.get(1)*2*2*16 + calo_geomid.get(2)*2*16 + calo_geomid.get(3)*16 + calo_geomid.get(4);
      else if (calo_geomid.get_type() == 1252) // gveto
	an_om_num = 520 + 128 + calo_geomid.get(1)*2*16 + calo_geomid.get(2)*16 + calo_geomid.get(3);

      if (om_num_to_kill[an_om_num])	
	CD.calorimeter_hits().erase(calo_hit);
      else ++calo_hit;
    }

  return PROCESS_SUCCESS;
}

void FLMANU_Killer::finalize()
{
  this->_set_initialized(false);
}
