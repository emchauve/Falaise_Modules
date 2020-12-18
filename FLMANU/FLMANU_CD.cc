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

uint32_t geomid_to_omnum(const geomtools::geom_id &geomid)
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

  // std::vector<unsigned short>   flag;
  std::vector<unsigned short> om_num;

  std::vector<float>          energy;
  std::vector<float>    sigma_energy;

  std::vector<float>            time;
  std::vector<float>      sigma_time;
  
  std::vector<float>     energy_true;
  std::vector<float>       time_true;
  std::vector<short>       nhit_true;
 
  std::string output_filename;
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
  
  if (myConfig.has_key("output_filename"))
    {
      output_filename = myConfig.fetch_string("output_filename");
    }

  ////////////////////////////////
  
  total_entries = 0; 
  selected_entries = 0;
  
  // to allow branching vectors in tree
  // gROOT->ProcessLine("#include<vector>");
  
  std::cout << "FLMANU_CD::initialize : output_filename = " << output_filename << std::endl;

  output_tree = new TTree ("output", "");
  output_tree->SetDirectory(0);

  // output_tree->Branch("flag",                  &flag);
  output_tree->Branch("om_num",              &om_num);

  output_tree->Branch("energy",              &energy);
  output_tree->Branch("sigma_energy",  &sigma_energy);
  output_tree->Branch("time",                  &time);
  output_tree->Branch("sigma_time",      &sigma_time);

  output_tree->Branch("energy_true",         &energy_true);
  output_tree->Branch("time_true",             &time_true);
  output_tree->Branch("nhit_true",             &nhit_true);

  this->_set_initialized(true);
}

dpp::chain_module::process_status FLMANU_CD::process(datatools::things &event)
{
  // std::cout << "FLMANU_CD::process()" << std::endl;

  int nhit_sd = 0;
  int nhit_cd = 0;

  ++total_entries;

  ////////
  // SD //
  ////////

  float tmp_energy_true[712];
  float tmp_time_true[712];
  short tmp_nhit_true[712];
  memset(tmp_energy_true, 0, 712*sizeof(float));
  memset(tmp_time_true, 0, 712*sizeof(float));
  memset(tmp_nhit_true, 0, 712*sizeof(short));

  mctools::simulated_data & SD = event.grab<mctools::simulated_data>("SD");

  std::vector<std::string> hit_categories;
  hit_categories.push_back("calo");
  hit_categories.push_back("xcalo");
  hit_categories.push_back("gveto");

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
    if (tmp_energy_true[om] > 0) nhit_sd++;

  if (nhit_sd <= 0)
    return dpp::base_module::PROCESS_STOP;

  ////////
  // CD //
  ////////

  float tmp_energy[712];
  float tmp_sigma_energy[712];
  float tmp_time[712];
  float tmp_sigma_time[712];
  memset(tmp_energy,       0, 712*sizeof(float));
  memset(tmp_sigma_energy, 0, 712*sizeof(float));
  memset(tmp_time,         0, 712*sizeof(float));
  memset(tmp_sigma_time,   0, 712*sizeof(float));

  const snemo::datamodel::calibrated_data & CD = event.get<snemo::datamodel::calibrated_data>("CD");

  for (const auto & calo_hit : CD.calibrated_calorimeter_hits())
    {
      short an_om_num = geomid_to_omnum(calo_hit->get_geom_id());

      if (an_om_num >= 712) {
	printf("*** om_num (cd) = %d\n", an_om_num);
	continue;}

      tmp_energy[an_om_num] = calo_hit->get_energy();
      tmp_sigma_energy[an_om_num] = calo_hit->get_sigma_energy();

      tmp_time[an_om_num] = calo_hit->get_time();
      tmp_sigma_time[an_om_num] = calo_hit->get_sigma_time();
    }
  
  for (int om=0; om<712; om++)
    if (tmp_energy[om] > 0) nhit_cd++;
  
  //

  // flag.clear();
  om_num.clear();

  energy.clear();
  sigma_energy.clear();
  time.clear();
  sigma_time.clear();

  energy_true.clear();
  time_true.clear();
  nhit_true.clear();

  for (int om=0; om<712; om++)
    {
      if (tmp_energy_true[om] <= 0)
	continue;

      om_num.push_back(om);

      energy.push_back(tmp_energy[om]);
      sigma_energy.push_back(tmp_sigma_energy[om]);
      time.push_back(tmp_time[om]);
      sigma_time.push_back(tmp_sigma_time[om]);

      energy_true.push_back(tmp_energy_true[om]);
      time_true.push_back(tmp_time_true[om]);
      nhit_true.push_back(tmp_nhit_true[om]);
    }
  
  output_tree->Fill();
  ++selected_entries;

  return dpp::base_module::PROCESS_SUCCESS;
}


void FLMANU_CD::finalize()
{
  std::cout << "FLMANU_CD::finalize   : " << selected_entries << "/" << total_entries << " selected" << std::endl;
  
  TFile *output_file = new TFile(output_filename.data(), "RECREATE");
  output_file->cd(); output_tree->Write(); output_file->Close();
}
