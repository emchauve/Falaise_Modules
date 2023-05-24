#include <bayeux/dpp/chain_module.h>
#include <bayeux/mctools/base_step_hit.h>
#include <bayeux/mctools/simulated_data.h>

#include <falaise/snemo/datamodels/event_header.h>
#include <falaise/snemo/datamodels/calibrated_data.h>
#include <falaise/snemo/datamodels/gg_track_utils.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <string>
#include <vector>
#include <array>

uint32_t geomid2omnum (const geomtools::geom_id &geomid)
{
  switch (geomid.get_type())
    {
    case 1302: // mwall
      return geomid.get(1)*20*13 + geomid.get(2)*13 + geomid.get(3);
    case 1232: // xwall
      return 520 + geomid.get(1)*2*2*16 + geomid.get(2)*2*16 + geomid.get(3)*16 + geomid.get(4);
    case 1252: // gveto
      return 520 + 128 + geomid.get(1)*2*16 + geomid.get(2)*16 + geomid.get(3);
    default:
      return -1;
    }
}


uint32_t geomid2cellnum (const geomtools::geom_id &geomid)
{
  switch (geomid.get_type())
    {
    case 1204: // gg
      return geomid.get(1)*1017 + geomid.get(3)*9 + geomid.get(2);
    default:
      return -1;
    }
}
  
//               geomid --> om_num
//
// mwall  min  0  0  0  0  =   0
// mwall  max  1  1 19 12  = 519
// 
// xwall  min  2  0 0 0  0 = 520
// xwall  max  3  1 1 1 15 = 647
//
// gveto  min  4  0 0  0   = 648
// gveto  min  5  1 1 15   = 711

class FLMANU_CD : public dpp::chain_module
{
public:
  FLMANU_CD();
  virtual ~FLMANU_CD();
    
  virtual void initialize (const datatools::properties &,
                           datatools::service_manager &,
			   dpp::module_handle_dict_type &);

  dpp::chain_module::process_status process (datatools::things & event);

  virtual void finalize ();
  
private:
  DPP_MODULE_REGISTRATION_INTERFACE(FLMANU_CD);

  unsigned long long total_entries;
  unsigned long long selected_entries;

  // unsigned int run_number;
  // unsigned int event_number;
  
  // enum {
  //   mwall_it_flag=0x01,
  //   mwall_fr_flag=0x02,
  //   xwall_it_flag=0x04,
  //   xwall_fr_flag=0x08,
  //   vwall_it_flag=0x10,
  //   vwall_fr_flag=0x20};

  unsigned long long event_id;

  std::vector<unsigned short> calo_om_num;
  std::vector<unsigned short> calo_flag;
  std::vector<float> calo_energy;
  std::vector<float> calo_energy_u;
  std::vector<float> calo_energy_bc;
  std::vector<float> calo_energy_bcu;
  std::vector<double> calo_time;
  //
  std::vector<float> calo_energy_true;
  std::vector<float> calo_energy_u_true;
  std::vector<float> calo_energy_bc_true;
  std::vector<float> calo_energy_bcu_true;
  std::vector<double> calo_time_true;

  //

  std::vector<unsigned short> tracker_cell_num;
  std::vector<unsigned short> tracker_flag;
  std::vector<float> tracker_anode_time;
  // std::vector<float> tracker_bottom_cathode_time;
  // std::vector<float> tracker_top_cathode_time;
  std::vector<float> tracker_r;
  std::vector<float> tracker_z;
  std::vector<float> tracker_drift_time_true; // in relative to avalanche beginning
  std::vector<float> tracker_anode_time_true; // in absolute
  std::vector<float> tracker_r_true;
  std::vector<float> tracker_z_true;

  //
 
  bool mw_only;
  bool mw8_only;

  float energy_cut;

  int nhit_calo_min;
  int nhit_tracker_min;

  std::string output_filename;

  TFile *output_file;
  TTree *output_tree;
};

DPP_MODULE_REGISTRATION_IMPLEMENT(FLMANU_CD, "FLMANU_CD")

FLMANU_CD::FLMANU_CD()
{
  // std::cout << "FLMANU_CD::FLMANU_CD()" << std::endl;
  output_filename = "flmanu-cd-output.root";
  this->_set_initialized(false);
}

FLMANU_CD::~FLMANU_CD()
{
  // std::cout << "FLMANU_CD::~FLMANU_CD" << std::endl;
  if (this->is_initialized()) this->finalize();
}

void FLMANU_CD::initialize (const datatools::properties & myConfig, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  // std::cout << "FLMANU_CD::initialize()" << std::endl;
  
  if (myConfig.has_key("mw_only"))
    mw_only = myConfig.fetch_boolean("mw_only");
  else mw_only = false;

  if (myConfig.has_key("mw8_only"))
    mw8_only = myConfig.fetch_boolean("mw8_only");
  else mw8_only = false;

  if (myConfig.has_key("energy_cut"))
    energy_cut = myConfig.fetch_real("energy_cut");
  else energy_cut = 0;

  if (myConfig.has_key("nhit_calo_min"))
    nhit_calo_min = myConfig.fetch_integer("nhit_calo_min");
  else nhit_calo_min = 0;

  if (myConfig.has_key("nhit_tracker_min"))
    nhit_tracker_min = myConfig.fetch_integer("nhit_tracker_min");
  else nhit_tracker_min = 0;

  if (myConfig.has_key("output_filename"))
    output_filename = myConfig.fetch_string("output_filename");

  ////////////////////////////////
  
  total_entries = 0;
  selected_entries = 0;
  
  // to allow branching vectors in tree
  // gROOT->ProcessLine("#include<vector>");
  
  std::cout << "FLMANU_CD::initialize : output_filename = " << output_filename << std::endl;

  output_file = new TFile(output_filename.data(), "RECREATE");

  output_tree = new TTree ("output", "");
  output_tree->SetDirectory(0);

  output_tree->Branch("event_id", &event_id);

  output_tree->Branch("calo.om_num",     &calo_om_num);
  // output_tree->Branch("calo.flag",       &calo_flag);
  output_tree->Branch("calo.energy",     &calo_energy);
  output_tree->Branch("calo.energy_u",   &calo_energy_u);
  output_tree->Branch("calo.energy_bc",  &calo_energy_bc);
  output_tree->Branch("calo.energy_bcu", &calo_energy_bcu);
  output_tree->Branch("calo.time",       &calo_time);
  //
  output_tree->Branch("calo.energy_true",     &calo_energy_true);
  output_tree->Branch("calo.energy_u_true",   &calo_energy_u_true);
  output_tree->Branch("calo.energy_bc_true",  &calo_energy_bc_true);
  output_tree->Branch("calo.energy_bcu_true", &calo_energy_bcu_true);
  output_tree->Branch("calo.time_true",       &calo_time_true);

  output_tree->Branch("tracker.cell_num",     &tracker_cell_num);
  // output_tree->Branch("tracker.flag",         &tracker_flag);
  output_tree->Branch("tracker.anode_time",   &tracker_anode_time);
  // output_tree->Branch("tracker.bottom_cathode_time",  &tracker_bottom_cathode_time);
  // output_tree->Branch("tracker.top_cathode_time",     &tracker_top_cathode_time);
  output_tree->Branch("tracker.r",            &tracker_r);
  output_tree->Branch("tracker.z",            &tracker_z);
  //
  output_tree->Branch("tracker.r_true",          &tracker_r_true);
  output_tree->Branch("tracker.z_true",          &tracker_z_true);
  output_tree->Branch("tracker.drift_time_true", &tracker_drift_time_true);
  output_tree->Branch("tracker.anode_time_true", &tracker_anode_time_true);

  this->_set_initialized(true);
}

dpp::chain_module::process_status FLMANU_CD::process(datatools::things &event)
{
  // std::cout << "FLMANU_CD::process()" << std::endl;

  int nhit_calo_sd = 0;
  int nhit_calo_cd = 0;

  int nhit_tracker_sd = 0;
  int nhit_tracker_cd = 0;

  event_id = total_entries++;

  ////////
  // SD //
  ////////

  mctools::simulated_data & SD = event.grab<mctools::simulated_data>("SD");

  // simulated calorimeter hit

  float tmp_calo_energy_true[712];
  float tmp_calo_energy_u_true[712];
  float tmp_calo_energy_bc_true[712];
  float tmp_calo_energy_bcu_true[712];
  double tmp_calo_time_true[712];

  memset(tmp_calo_energy_true, 0, 712*sizeof(float));
  memset(tmp_calo_energy_u_true, 0, 712*sizeof(float));
  memset(tmp_calo_energy_bc_true, 0, 712*sizeof(float));
  memset(tmp_calo_energy_bcu_true, 0, 712*sizeof(float));
  memset(tmp_calo_time_true, 0, 712*sizeof(double));

  std::vector<std::string> hit_categories;
  hit_categories.push_back("calo");

  if (!(mw_only || mw8_only))
    {
      hit_categories.push_back("xcalo");
      hit_categories.push_back("gveto");
    }
  
  for (std::string hit_category : hit_categories)
    {
      if (!SD.has_step_hits(hit_category))
	continue;

      for (auto & a_step_hit : SD.get_step_hits(hit_category))
	{
	  unsigned short an_om_num = geomid2omnum(a_step_hit->get_geom_id());

	  if (an_om_num >= 712) {
	    printf("*** om_num (sd) = %d\n", an_om_num);
	    continue;}

	  if (tmp_calo_energy_true[an_om_num] == 0) // first hit in this OM
	    tmp_calo_time_true[an_om_num] = a_step_hit->get_time_start();

	  else if (a_step_hit->get_time_start() < tmp_calo_time_true[an_om_num])
	    tmp_calo_time_true[an_om_num] = a_step_hit->get_time_start();
	  
	  tmp_calo_energy_true[an_om_num] += a_step_hit->get_energy_deposit();
	}
    }
  
  for (int om=0; om<712; om++)
    if (tmp_calo_energy_true[om] > 0) nhit_calo_sd++;

  // if (nhit_calo_sd <= 0)
  //   return dpp::base_module::PROCESS_STOP;

  // simulated tracker hit

  double tmp_tracker_r_true[2034];
  double tmp_tracker_z_true[2034];
  double tmp_tracker_drift_time_true[2034];
  double tmp_tracker_anode_time_true[2034];
  memset(tmp_tracker_r_true,    0, 2034*sizeof(double));
  memset(tmp_tracker_z_true,    0, 2034*sizeof(double));
  memset(tmp_tracker_drift_time_true, 0, 2034*sizeof(double));
  memset(tmp_tracker_anode_time_true, 0, 2034*sizeof(double));

  if (SD.has_step_hits("gg"))
    {
      for (auto & a_step_hit : SD.get_step_hits("gg"))
	{
	  unsigned short a_cell_num = geomid2cellnum(a_step_hit->get_geom_id());

	  const double z_true = a_step_hit->get_position_start().z();

	  // do not trigger the cell if track is inside cell
	  if (std::abs(z_true) > 1399 * CLHEP::millimeter)
	    continue;

	  const double step_hit_time = a_step_hit->get_time_start();

	  const datatools::properties & step_hit_aux = a_step_hit->get_auxiliaries();
	  const std::string & r_true_aux_key = snemo::datamodel::gg_track::minimum_approach_distance_key();

	  const double r_true = step_hit_aux.fetch_real(r_true_aux_key);
	  const double r_true_norm = r_true / (44.0 * CLHEP::millimeter);

	  const double drift_time_1 = (std::exp(r_true_norm/2.44549e-01)-1)/1.69542e+00;
	  const double drift_time_2 = (std::exp(r_true_norm/5.35294e-02)-1)/3.13692e+03;
	  const double drift_time = std::max(drift_time_1, drift_time_2) * CLHEP::microsecond;

	  const double drift_time_true = drift_time;
	  const double anode_time_true = step_hit_time + drift_time;

	  if (tmp_tracker_drift_time_true[a_cell_num] == 0)
	    {
	      // first step hit in this cell
	      tmp_tracker_r_true[a_cell_num] = r_true;
	      tmp_tracker_z_true[a_cell_num] = z_true;
	      tmp_tracker_drift_time_true[a_cell_num] = drift_time_true;
	      tmp_tracker_anode_time_true[a_cell_num] = anode_time_true;
	    }
	  else if (anode_time_true < tmp_tracker_anode_time_true[a_cell_num])
	    {
	      // another step hit in this cell: replace true information
	      // check if the anode start time is ealier (taking into account
	      // step hit time + anode drift time = "anode time")

	      // printf("+++ updating drift time %7.1f ns => %7.1f (cell %d in event %d)\n",
	      // 	     tmp_tracker_anode_time_true[a_cell_num], drift_time_true, a_cell_num, event_id);

	      tmp_tracker_r_true[a_cell_num] = r_true;
	      tmp_tracker_z_true[a_cell_num] = z_true;
	      tmp_tracker_drift_time_true[a_cell_num] = drift_time_true;
	      tmp_tracker_anode_time_true[a_cell_num] = anode_time_true;
	    }
	}
    }

  ////////
  // CD //
  ////////

  const snemo::datamodel::calibrated_data & CD = event.get<snemo::datamodel::calibrated_data>("CD");

  // process calibrated calorimeter hits

  float tmp_calo_energy[712];
  float tmp_calo_energy_u[712];
  float tmp_calo_energy_bc[712];
  float tmp_calo_energy_bcu[712];
  double tmp_calo_time[712];

  memset(tmp_calo_energy,     0, 712*sizeof(float));
  memset(tmp_calo_energy_u,   0, 712*sizeof(float));
  memset(tmp_calo_energy_bc,  0, 712*sizeof(float));
  memset(tmp_calo_energy_bcu, 0, 712*sizeof(float));
  memset(tmp_calo_time,       0, 712*sizeof(double));

  for (const auto & calo_hit : CD.calorimeter_hits())
    {
      const unsigned short an_om_num = geomid2omnum(calo_hit->get_geom_id());

      const datatools::properties& calo_hit_aux = calo_hit->get_auxiliaries();

      if (calo_hit_aux.has_key("edep"))
	{
	  if (tmp_calo_energy[an_om_num] > 0)
	    printf("*** calibrated calorimeter hit already set for om %d in event_id %d\n", an_om_num, event_id);

	  // tmp_calo_energy_true[an_om_num] = calo_hit_aux.fetch_real("edep");
	  tmp_calo_energy_u_true[an_om_num] = calo_hit_aux.fetch_real("edep_u");
	  tmp_calo_energy_bc_true[an_om_num] = calo_hit_aux.fetch_real("edep_bc");
	  tmp_calo_energy_bcu_true[an_om_num] = calo_hit_aux.fetch_real("edep_bcu");

	  tmp_calo_energy[an_om_num] = calo_hit_aux.fetch_real("evis");
	  tmp_calo_energy_u[an_om_num] = calo_hit_aux.fetch_real("evis_u");
	  tmp_calo_energy_bc[an_om_num] = calo_hit_aux.fetch_real("evis_bc");
	  tmp_calo_energy_bcu[an_om_num] = calo_hit_aux.fetch_real("evis_bcu");
	}

      tmp_calo_time[an_om_num] = calo_hit->get_time();

      nhit_calo_cd++;
    }

  // process calibrated tracker hits

  float tmp_tracker_anode_time[2034];
  // float tmp_tracker_bottom_cathode_time[2034];
  // float tmp_tracker_top_cathode_time[2034];
  float tmp_tracker_r[2034];
  float tmp_tracker_z[2034];

  memset(tmp_tracker_anode_time,   0, 2034*sizeof(float));
  // memset(tmp_tracker_bottom_cathode_time, 0, 2034*sizeof(float));
  // memset(tmp_tracker_top_cathode_time,    0, 2034*sizeof(float));
  memset(tmp_tracker_r,            0, 2034*sizeof(float));
  memset(tmp_tracker_z,            0, 2034*sizeof(float));

  for (const auto & tracker_hit : CD.tracker_hits())
    {
      // short a_cell_num = geomid2cellnum(hit->get_geom_id());

      unsigned short a_cell_num = 9 * 113 * tracker_hit->get_side();
      a_cell_num += 9 * tracker_hit->get_row();
      a_cell_num += tracker_hit->get_layer();

      if (tmp_tracker_r[a_cell_num] > 0)
	    printf("*** calibrated tracker hit already set for cell %d in event_id %d\n", a_cell_num, event_id);

      tmp_tracker_anode_time[a_cell_num] = tracker_hit->get_anode_time();
      tmp_tracker_r[a_cell_num] = tracker_hit->get_r();
      tmp_tracker_z[a_cell_num] = tracker_hit->get_z();
    }

  // reset & fill calo data for output tree

  calo_om_num.clear();
  calo_flag.clear();
  calo_energy.clear();
  calo_energy_u.clear();
  calo_energy_bc.clear();
  calo_energy_bcu.clear();
  calo_time.clear();

  calo_energy_true.clear();
  calo_energy_u_true.clear();
  calo_energy_bc_true.clear();
  calo_energy_bcu_true.clear();
  calo_time_true.clear();

  for (unsigned short om=0; om<712; om++)
    {
      // if (tmp_calo_energy_true[om] <= 0)
      if (tmp_calo_energy[om] <= 0)
	continue;

      if (tmp_calo_energy[om] < energy_cut)
	continue;

      if (mw8_only && ((om%13)==0 || (om%13)==12))
	continue;

      calo_om_num.push_back(om);
      calo_energy.push_back(tmp_calo_energy[om]);
      calo_energy_u.push_back(tmp_calo_energy_u[om]);
      calo_energy_bc.push_back(tmp_calo_energy_bc[om]);
      calo_energy_bcu.push_back(tmp_calo_energy_bcu[om]);
      calo_time.push_back(tmp_calo_time[om]);

      calo_energy_true.push_back(tmp_calo_energy_true[om]);
      calo_energy_u_true.push_back(tmp_calo_energy_u_true[om]);
      calo_energy_bc_true.push_back(tmp_calo_energy_bc_true[om]);
      calo_energy_bcu_true.push_back(tmp_calo_energy_bcu_true[om]);
      calo_time_true.push_back(tmp_calo_time_true[om]);
    }

  // reset & fill tracker data for output tree

  tracker_cell_num.clear();
  tracker_flag.clear();
  tracker_anode_time.clear();
  // tracker_bottom_cathode_time.clear();
  // tracker_top_cathode_time.clear();
  tracker_r.clear();
  tracker_z.clear();

  tracker_drift_time_true.clear();
  tracker_anode_time_true.clear();
  tracker_r_true.clear();
  tracker_z_true.clear();

  for (unsigned short cell=0; cell<2034; cell++)
    {
      if (tmp_tracker_r_true[cell] <= 0)
	// if (tmp_tracker_r[cell] <= 0)
	continue;

      tracker_cell_num.push_back(cell);
      tracker_anode_time.push_back(tmp_tracker_anode_time[cell]);
      tracker_r.push_back(tmp_tracker_r[cell]);
      tracker_z.push_back(tmp_tracker_z[cell]);

      tracker_drift_time_true.push_back(tmp_tracker_drift_time_true[cell]);
      tracker_anode_time_true.push_back(tmp_tracker_anode_time_true[cell]);
      tracker_r_true.push_back(tmp_tracker_r_true[cell]);
      tracker_z_true.push_back(tmp_tracker_z_true[cell]);

      nhit_tracker_cd++;
    }

  /////////////////////////////

  if ((nhit_calo_cd >= nhit_calo_min) && (nhit_tracker_cd >= nhit_tracker_min))
    {
      output_tree->Fill();
      ++selected_entries;
    }

  return dpp::base_module::PROCESS_SUCCESS;
}


void FLMANU_CD::finalize()
{
  std::cout << "FLMANU_CD::finalize   : " << selected_entries << "/" << total_entries << " selected" << std::endl;
  
  output_file->cd();
  output_tree->Write();
  output_file->Close();
}
