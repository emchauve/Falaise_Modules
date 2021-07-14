#include <bayeux/dpp/chain_module.h>
#include <bayeux/mctools/base_step_hit.h>
#include <bayeux/mctools/simulated_data.h>

#include <falaise/snemo/datamodels/event_header.h>
#include <falaise/snemo/datamodels/calibrated_data.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <string>
#include <vector>
#include <array>

uint32_t geomid_to_omnum (const geomtools::geom_id &geomid)
{
  switch (geomid.get_type()) {
  case 1302: // mwall
    return geomid.get(1)*20*13 + geomid.get(2)*13 + geomid.get(3);
  case 1232: // xwall
    return 520 + geomid.get(1)*2*2*16 + geomid.get(2)*2*16 + geomid.get(3)*16 + geomid.get(4);
  case 1252: // gveto
    return 520 + 128 + geomid.get(1)*2*16 + geomid.get(2)*16 + geomid.get(3);
  default: return -1;}
}


// uint32_t geomid_to_cellnum (const geomtools::geom_id &geomid)
// {
//   switch (geomid.get_type()) {
//   case 1302: // mwall
//     return geomid.get(1)*20*13 + geomid.get(2)*13 + geomid.get(3);
//   case 1232: // xwall
//     return 520 + geomid.get(1)*2*2*16 + geomid.get(2)*2*16 + geomid.get(3)*16 + geomid.get(4);
//   case 1252: // gveto
//     return 520 + 128 + geomid.get(1)*2*16 + geomid.get(2)*16 + geomid.get(3);
//   default: return -1;}
// }
  
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

  // std::vector<unsigned short>   flag;
  std::vector<unsigned short> om_num;

  std::vector<float> energy;
  std::vector<float> energy_u;
  std::vector<float> energy_bc;
  std::vector<float> energy_bcu;
  std::vector<float> sigma_energy;
  std::vector<float> sigma_energy_u;
  std::vector<float> sigma_energy_bc;
  std::vector<float> sigma_energy_bcu;
  std::vector<float> time;
  std::vector<float> sigma_time;
  
  std::vector<float> energy_true;
  std::vector<float> energy_u_true;
  std::vector<float> energy_bc_true;
  std::vector<float> energy_bcu_true;
  std::vector<float> time_true;
  std::vector<short> nhit_true;

  //

  std::vector<unsigned short> tracker_cell_num;

  std::vector<float> tracker_r;
  std::vector<float> tracker_z;
  std::vector<float> tracker_anode_time;
  std::vector<float> tracker_delayed_time;

  //
 
  bool mw_only;
  bool mw8_only;

  float energy_cut;

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

  // output_tree->Branch("flag",                  &flag);
  output_tree->Branch("om_num", &om_num);

  output_tree->Branch("energy", &energy);
  output_tree->Branch("energy_u", &energy_u);
  output_tree->Branch("energy_bc", &energy_bc);
  output_tree->Branch("energy_bcu", &energy_bcu);

  // output_tree->Branch("sigma_energy", &sigma_energy);
  // output_tree->Branch("sigma_energy_u", &sigma_energy_u);
  // output_tree->Branch("sigma_energy_bc", &sigma_energy_bc);
  // output_tree->Branch("sigma_energy_bcu", &sigma_energy_bcu);

  output_tree->Branch("time",                  &time);
  // output_tree->Branch("sigma_time",      &sigma_time);

  output_tree->Branch("energy_true",  &energy_true);
  output_tree->Branch("energy_u_true",  &energy_u_true);
  output_tree->Branch("energy_bc_true",  &energy_bc_true);
  output_tree->Branch("energy_bcu_true",  &energy_bcu_true);
  output_tree->Branch("time_true",      &time_true);
  output_tree->Branch("nhit_true",      &nhit_true);


  output_tree->Branch("tracker_cell_num",     &tracker_cell_num);
  output_tree->Branch("tracker_r",            &tracker_r);
  output_tree->Branch("tracker_z",            &tracker_z);
  output_tree->Branch("tracker_anode_time",   &tracker_anode_time);
  output_tree->Branch("tracker_delayed_time", &tracker_delayed_time);

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

  float tmp_energy_true[712];
  float tmp_energy_u_true[712];
  float tmp_energy_bc_true[712];
  float tmp_energy_bcu_true[712];
  float tmp_time_true[712];
  short tmp_nhit_true[712];

  memset(tmp_energy_true, 0, 712*sizeof(float));
  memset(tmp_energy_u_true, 0, 712*sizeof(float));
  memset(tmp_energy_bc_true, 0, 712*sizeof(float));
  memset(tmp_energy_bcu_true, 0, 712*sizeof(float));
  memset(tmp_time_true, 0, 712*sizeof(float));
  memset(tmp_nhit_true, 0, 712*sizeof(short));

  mctools::simulated_data & SD = event.grab<mctools::simulated_data>("SD");

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
	  short an_om_num = geomid_to_omnum(a_step_hit->get_geom_id());

	  if (an_om_num >= 712) {
	    printf("*** om_num (sd) = %d\n", an_om_num);
	    continue;}

	  if (tmp_energy_true[an_om_num] == 0) // first hit in this OM
	    tmp_time_true[an_om_num] = a_step_hit->get_time_start();

	  else if (a_step_hit->get_time_start() < tmp_time_true[an_om_num])
	    tmp_time_true[an_om_num] = a_step_hit->get_time_start();
	  
	  tmp_energy_true[an_om_num] += a_step_hit->get_energy_deposit();
	  
	  tmp_nhit_true[an_om_num];
	}
    }
  
  for (int om=0; om<712; om++)
    if (tmp_energy_true[om] > 0) nhit_calo_sd++;

  // if (nhit_calo_sd <= 0)
  //   return dpp::base_module::PROCESS_STOP;


  ////////
  // CD //
  ////////

  float tmp_energy[712];
  float tmp_energy_u[712];
  float tmp_energy_bc[712];
  float tmp_energy_bcu[712];
  float tmp_sigma_energy[712];
  float tmp_time[712];
  float tmp_sigma_time[712];

  memset(tmp_energy, 0, 712*sizeof(float));
  memset(tmp_energy_u, 0, 712*sizeof(float));
  memset(tmp_energy_bc, 0, 712*sizeof(float));
  memset(tmp_energy_bcu, 0, 712*sizeof(float));
  memset(tmp_sigma_energy, 0, 712*sizeof(float));
  memset(tmp_time,         0, 712*sizeof(float));
  memset(tmp_sigma_time,   0, 712*sizeof(float));

  const snemo::datamodel::calibrated_data & CD = event.get<snemo::datamodel::calibrated_data>("CD");

  for (const auto & calo_hit : CD.calorimeter_hits())
    {
      short an_om_num = geomid_to_omnum(calo_hit->get_geom_id());

      if (an_om_num >= 712) {
	printf("*** om_num (cd) = %d\n", an_om_num);
	continue;}

      tmp_energy[an_om_num] = calo_hit->get_energy();
      tmp_sigma_energy[an_om_num] = calo_hit->get_sigma_energy();

      // retrieve energy detail in auxilliaries

      const datatools::properties& calo_hit_aux = calo_hit->get_auxiliaries();

      if (calo_hit_aux.has_key("edep"))
	{
	  // tmp_energy_true[an_om_num] = calo_hit_aux.fetch_real("edep");
	  tmp_energy_u_true[an_om_num] = calo_hit_aux.fetch_real("edep_u");
	  tmp_energy_bc_true[an_om_num] = calo_hit_aux.fetch_real("edep_bc");
	  tmp_energy_bcu_true[an_om_num] = calo_hit_aux.fetch_real("edep_bcu");

	  // tmp_energy[an_om_num] = calo_hit_aux.fetch_real("evis");
	  tmp_energy_u[an_om_num] = calo_hit_aux.fetch_real("evis_u");
	  tmp_energy_bc[an_om_num] = calo_hit_aux.fetch_real("evis_bc");
	  tmp_energy_bcu[an_om_num] = calo_hit_aux.fetch_real("evis_bcu");
	}

      tmp_time[an_om_num] = calo_hit->get_time();
      tmp_sigma_time[an_om_num] = calo_hit->get_sigma_time();
    }
  
  for (int om=0; om<712; om++)
    if (tmp_energy[om] > 0) nhit_calo_cd++;
  
  // flag.clear();
  om_num.clear();

  energy.clear();
  energy_u.clear();
  energy_bc.clear();
  energy_bcu.clear();
  sigma_energy.clear();
  sigma_energy_u.clear();
  sigma_energy_bc.clear();
  sigma_energy_bcu.clear();

  time.clear();
  sigma_time.clear();

  energy_true.clear();
  energy_u_true.clear();
  energy_bc_true.clear();
  energy_bcu_true.clear();
  time_true.clear();
  nhit_true.clear();

  for (unsigned short om=0; om<712; om++)
    {
      if (tmp_energy_true[om] <= 0)
	continue;

      if (tmp_energy[om] < energy_cut)
	continue;

      if (mw8_only && ((om%13)==0 || (om%13)==12))
	continue;

      om_num.push_back(om);

      energy.push_back(tmp_energy[om]);
      energy_u.push_back(tmp_energy_u[om]);
      energy_bc.push_back(tmp_energy_bc[om]);
      energy_bcu.push_back(tmp_energy_bcu[om]);
      sigma_energy.push_back(tmp_sigma_energy[om]);
      time.push_back(tmp_time[om]);
      sigma_time.push_back(tmp_sigma_time[om]);

      energy_true.push_back(tmp_energy_true[om]);
      energy_u_true.push_back(tmp_energy_u_true[om]);
      energy_bc_true.push_back(tmp_energy_bc_true[om]);
      energy_bcu_true.push_back(tmp_energy_bcu_true[om]);
      time_true.push_back(tmp_time_true[om]);
      nhit_true.push_back(tmp_nhit_true[om]);
    }
  
  /////////////////////////////
  // calibrated tracker hits //
  /////////////////////////////

  float tmp_tracker_r[2034];
  float tmp_tracker_z[2034];
  float tmp_tracker_anode_time[2034];
  float tmp_tracker_delayed_time[2034];

  memset(tmp_tracker_r,            0, 2034*sizeof(float));
  memset(tmp_tracker_z,            0, 2034*sizeof(float));
  memset(tmp_tracker_anode_time,   0, 2034*sizeof(float));
  memset(tmp_tracker_delayed_time, 0, 2034*sizeof(float));

  for (const auto & tracker_hit : CD.tracker_hits())
    {
      // short a_cell_num = geomid_to_cellnum(hit->get_geom_id());

      unsigned short a_cell_num = 9 * 113 * tracker_hit->get_side();
      a_cell_num += 9 * tracker_hit->get_row();
      a_cell_num += tracker_hit->get_layer();

      if (a_cell_num >= 2034)
	printf("cell num = %4d for cell %d/%d/%03d\n", a_cell_num, tracker_hit->get_side(), tracker_hit->get_row(), tracker_hit->get_layer());

      tmp_tracker_r[a_cell_num] = tracker_hit->get_r();
      tmp_tracker_z[a_cell_num] = tracker_hit->get_z();
      tmp_tracker_anode_time[a_cell_num] = tracker_hit->get_anode_time();
      tmp_tracker_delayed_time[a_cell_num] = tracker_hit->get_delayed_time();
    }

  //

  tracker_cell_num.clear();

  tracker_r.clear();
  tracker_z.clear();
  tracker_anode_time.clear();
  tracker_delayed_time.clear();

  for (unsigned short cell=0; cell<2034; cell++)
    {
      if (tmp_tracker_r[cell] == 0)
	continue;

      tracker_cell_num.push_back(cell);

      tracker_r.push_back(tmp_tracker_r[cell]);
      tracker_z.push_back(tmp_tracker_z[cell]);
      tracker_anode_time.push_back(tmp_tracker_anode_time[cell]);
      tracker_delayed_time.push_back(tmp_tracker_delayed_time[cell]);

      nhit_tracker_cd++;
    }

  /////////////////////////////

  output_tree->Fill();
  ++selected_entries;

  return dpp::base_module::PROCESS_SUCCESS;
}


void FLMANU_CD::finalize()
{
  std::cout << "FLMANU_CD::finalize   : " << selected_entries << "/" << total_entries << " selected" << std::endl;
  
  output_file->cd();
  output_tree->Write();
  output_file->Close();
}
