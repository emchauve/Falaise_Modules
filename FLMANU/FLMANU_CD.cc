#include <bayeux/dpp/chain_module.h>
// #include <falaise/snemo/processing/module.h>

#include <falaise/snemo/datamodels/event_header.h>
#include <falaise/snemo/datamodels/calibrated_data.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <string>
#include <vector>
#include <array>

uint32_t geomid_to_pmt(const geomtools::geom_id &geomid) {
  switch (geomid.get_type()) {
  case 1302: // mwall
    return geomid.get(1)*20*13 + geomid.get(2)*13 + geomid.get(3);
  case 1232: // xwall
    return 520 + geomid.get(1)*2*2*16 + geomid.get(2)*2*16 + geomid.get(3)*16 + geomid.get(4);
  case 1252: // gveto
    return 520 + 128 + geomid.get(1)*2*16 + geomid.get(2)*16 + geomid.get(3);
  default: return -1;}
}
  
//               geomid --> pmt
//
// mwall  min  0 0-0 00  =  000
// mwall  max  1 1-9 12  =  519
// 
// xwall  min  2 0 0 00  =  520
// xwall  max  3 1 8 12  =  647
//
// gveto  min  4 0 0 00  =  648
// gveto  min  5 1 1 15  =  711


class FLMANU_CD : public dpp::chain_module
{
public:
  FLMANU_CD();
  virtual ~FLMANU_CD();
  
  // FLMANU_CD (const datatools::properties & myConfig,
  // 	    const datatools::service_manager & flServices);
  
  virtual void initialize (const datatools::properties &,
                           datatools::service_manager &,
			   dpp::module_handle_dict_type &);

  dpp::chain_module::process_status process (datatools::things & workItem);
  // falaise::processing::status process (datatools::things& workItem);

  virtual void finalize ();
  
private:
  DPP_MODULE_REGISTRATION_INTERFACE(FLMANU_CD);

  unsigned long long total_entries;
  unsigned long long selected_entries;

  // unsigned int run_number;
  // unsigned int event_number;
  
  enum {
    mwall_it_flag=0x01,
    mwall_fr_flag=0x02,
    xwall_it_flag=0x04,
    xwall_fr_flag=0x08,
    vwall_it_flag=0x10,
    vwall_fr_flag=0x20};

  unsigned short nhit;

  std::vector<unsigned short>   flag;
  std::vector<unsigned short> geomid;
  std::vector<float>          energy;
  std::vector<float>    sigma_energy;
  std::vector<float>            time;
  std::vector<float>      sigma_time;

  std::vector<float>     energy_true;
  std::vector<float>       time_true;
  
  // unsigned short        nhit0;
  // std::vector<float>  energy0;
  // std::vector<float>    time0;
 
  std::string output_filename;
  TTree *output_tree;

  unsigned int min_nhit;

};

DPP_MODULE_REGISTRATION_IMPLEMENT(FLMANU_CD, "FLMANU_CD")
// FALAISE_REGISTER_MODULE(FLMANU_CD);

FLMANU_CD::FLMANU_CD()
{
  //std::cout << "FLMANU_CD::FLMANU_CD()" << std::endl;
  output_filename = "output.root";
  this->_set_initialized(false);
}

FLMANU_CD::~FLMANU_CD()
{
  // std::cout << "FLMANU_CD::~FLMANU_CD" << std::endl;
  if (this->is_initialized()) this->finalize();
}

void FLMANU_CD::initialize (const datatools::properties & myConfig, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  // std::cout << "FLMANU_CD::initialize(...)" << std::endl;
  
  if (myConfig.has_key("output_filename"))
    {
      output_filename = myConfig.fetch_string("output_filename");
      std::cout << "FLMANU_CD::initialize : setting output_filename = " << output_filename << std::endl;
    }

  if (myConfig.has_key("min_nhit"))
    {
      min_nhit = myConfig.fetch_integer("min_nhit");
      std::cout << "FLMANU_CD::initialize : setting min_nhit = " << min_nhit << std::endl;
    }
  
  else min_nhit = 0;  

  ////////////////////////////////
  
  total_entries = 0; 
  selected_entries = 0;
  
  // to allow branching vectors in tree
  // gROOT->ProcessLine("#include<vector>");
  
  output_tree = new TTree ("output", "");
  output_tree->SetDirectory(0);

  // output_tree->Branch("run", &run_number);
  // output_tree->Branch("event", &event_number);

  output_tree->Branch("event",             &total_entries);

  output_tree->Branch("nhit", &nhit);

  output_tree->Branch("flag",                  &flag);
  output_tree->Branch("geomid",              &geomid);
  output_tree->Branch("energy",              &energy);
  output_tree->Branch("sigma_energy",  &sigma_energy);
  output_tree->Branch("time",                  &time);
  output_tree->Branch("sigma_time",      &sigma_time);

  output_tree->Branch("energy_true",         &energy_true);
  output_tree->Branch("time_true",             &time_true);

  // output_tree->Branch("nhit0",   &nhit0;
  // output_tree->Branch("energy0", &energy0);
  // output_tree->Branch("time0", &time0);

  this->_set_initialized(true);
}

// falaise::processing::status FLMANU_CD::process(datatools::things &workItem)
dpp::chain_module::process_status FLMANU_CD::process(datatools::things &workItem)
{
  // if ((++total_entries % 1000000) == 0)
  //   DT_LOG_NOTICE(get_logging_priority(), Form("[%09d] events", total_entries));
  
  nhit = 0;  
  flag.clear();
  geomid.clear();
  energy.clear();
  sigma_energy.clear();
  time.clear();
  sigma_time.clear();
  energy_true.clear();
  time_true.clear();

  ++total_entries;
  
  // nhit0 = 0;
  // energy0.clear();
  // time0.clear();
  
  const snemo::datamodel::calibrated_data & CD = workItem.get<snemo::datamodel::calibrated_data>("CD");

  if (CD.calibrated_calorimeter_hits().size() == 0)
    return PROCESS_STOP;
    // return falaise::processing::status::PROCESS_STOP;

  for (const auto & calo_hit : CD.calibrated_calorimeter_hits())
    {
      unsigned short _flag_ = 0;

      const geomtools::geom_id & calo_hit_geomid = calo_hit->get_geom_id();
	
      switch (calo_hit_geomid.get_type())
	{
	case 1302: // main wall
	  if (calo_hit_geomid.get(1) == 0)      _flag_ |= mwall_it_flag;
	  else if (calo_hit_geomid.get(1) == 1) _flag_ |= mwall_fr_flag;
	  else std::cout << "*** unexpected main wall OM type = " << calo_hit_geomid.get(0) << std::endl;
	  break;

	case 1232: // xwall
	  if (calo_hit_geomid.get(1) == 0)      _flag_ |= xwall_it_flag;
	  else if (calo_hit_geomid.get(1) == 1) _flag_ |= xwall_fr_flag;
	  else std::cout << "*** unexpected X-wall OM type = " << calo_hit_geomid.get(0) << std::endl; 
	  break;

	case 1252: // vwall
	  if (calo_hit_geomid.get(1) == 0)      _flag_ |= vwall_it_flag;
	  else if (calo_hit_geomid.get(1) == 1) _flag_ |= vwall_fr_flag;
	  else std::cout << "*** unexpected V-wall OM type = " << calo_hit_geomid.get(0) << std::endl;
	  break;
	
	default:
	  std::cout << "*** unexpected geom_id = " << calo_hit_geomid.get_type()<< std::endl;
	  break;
	}

      flag.push_back(_flag_);

      geomid.push_back(geomid_to_pmt(calo_hit_geomid));

      energy.push_back(calo_hit->get_energy());
      sigma_energy.push_back(calo_hit->get_sigma_energy()*2.35482);

      time.push_back(calo_hit->get_time());
      sigma_time.push_back(calo_hit->get_sigma_time());

      // energy_true.push_back();
      // time_true.push_back();

      ++nhit;
    }
  
  if (nhit < min_nhit)
    return PROCESS_STOP;
    // return falaise::processing::status::PROCESS_STOP;

  output_tree->Fill();
  ++selected_entries;

  return PROCESS_OK;
  // return falaise::processing::status::PROCESS_OK;
}


void FLMANU_CD::finalize()
{
  std::cout << "FLMANU_CD::finalize   : " << selected_entries << "/" << total_entries << " selected" << std::endl;

  TFile *output_file = new TFile(output_filename.data(), "RECREATE");
  output_file->cd(); output_tree->Write(); output_file->Close();
}
