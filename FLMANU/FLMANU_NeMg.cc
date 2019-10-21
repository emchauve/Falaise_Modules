#include <bayeux/dpp/chain_module.h>
#include <bayeux/genbb_help/primary_event.h>
#include <bayeux/mctools/base_step_hit.h>
#include <bayeux/mctools/simulated_data.h>

#include <falaise/snemo/datamodels/event_header.h>
#include <falaise/snemo/datamodels/calibrated_data.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/datamodels/particle_track.h>
#include <falaise/snemo/datamodels/tracker_trajectory.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include <falaise/snemo/datamodels/tracker_trajectory_solution.h>

#include <falaise/snemo/datamodels/helix_trajectory_pattern.h>
#include <falaise/snemo/datamodels/line_trajectory_pattern.h>

#include <gsl/gsl_cdf.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <string>
#include <vector>
#include <array>

namespace tof_tools
{
  const float kC  = 299.792458;   // mm/ns
  const float kMe = 0.5109989461; // MeV

  inline float beta (float energy, float mass) {
    if (mass == 0.0) return 1.0;
    else return (std::sqrt(energy*(energy+2.*mass)))/(energy+mass);}
  
  inline float tof_th (float energy, float trlen, float mass) {
    return trlen/(beta(energy,mass)*kC);}

  inline float sigma2_tof_th (float tof, float energy, float sigma_e, float mass) {
    return std::pow((tof*mass*mass)/(energy*(energy+mass)*(energy+2*mass)),2)*energy*sigma_e*sigma_e;}

};


class FLMANU_NeMg : public dpp::chain_module
{
public:
  FLMANU_NeMg();
  virtual ~FLMANU_NeMg();
  virtual void initialize(const datatools::properties &myConfig,
                          datatools::service_manager& flServices,
                          dpp::module_handle_dict_type& moduleDict);

  virtual dpp::chain_module::process_status process(datatools::things &event);
  virtual void finalize();

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

  
private:
  DPP_MODULE_REGISTRATION_INTERFACE(FLMANU_NeMg);

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

  unsigned short e;
  unsigned short g;

  std::vector<unsigned short>   flag_e;
  // std::vector<unsigned short> ncalo_e;
  std::vector<unsigned short> geomid_e;
  std::vector<float>          energy_e;
  std::vector<float>    sigma_energy_e;
  std::vector<float>            time_e;
  std::vector<float>      sigma_time_e;
  std::vector<float>     time_th_int_e;
  std::vector<float>     time_th_ext_e;
  std::vector<float>    track_length_e;

  std::vector<float>          x_foil_e;
  std::vector<float>          y_foil_e;
  std::vector<float>          z_foil_e;

  std::vector<float>          x_calo_e;
  std::vector<float>          y_calo_e;
  std::vector<float>          z_calo_e;

  std::vector<float>          energy_e_true;
  std::vector<float>            time_e_true;

  // std::vector<float>          x_foil_e_true;
  // std::vector<float>          y_foil_e_true;
  // std::vector<float>          z_foil_e_true;

  std::vector<float>          x_calo_e_true;
  std::vector<float>          y_calo_e_true;
  std::vector<float>          z_calo_e_true;

  std::vector<unsigned short>   flag_g;
  // std::vector<unsigned short> ncalo_g;
  std::vector<unsigned short> geomid_g;
  std::vector<float>          energy_g;
  std::vector<float>    sigma_energy_g;
  std::vector<float>            time_g;
  std::vector<float>      sigma_time_g;

  std::vector<float>          energy_g_true;
  std::vector<float>            time_g_true;

  unsigned short e0;
  unsigned short g0;

  std::vector<float>        energy_e0_true;
  std::vector<float>        energy_g0_true;

  float                            x0_true;
  float                            y0_true;
  float                            z0_true;

  float pint_e1e2;
  float pext_e1e2;

  float pint_e1g1;
  float pext_e1g1;

  float pint_e2g1;
  float pext_e2g1;



  TTree *output_tree;

  std::string output_filename;

  unsigned int min_e;
  unsigned int min_g;

};

DPP_MODULE_REGISTRATION_IMPLEMENT(FLMANU_NeMg, "FLMANU_NeMg");

FLMANU_NeMg::FLMANU_NeMg() : dpp::chain_module()
{
  output_filename = "output.root";
}

FLMANU_NeMg::~FLMANU_NeMg()
{
  this->finalize();
}

void FLMANU_NeMg::initialize(const datatools::properties  &myConfig,
			      datatools::service_manager   &/*&flServices*/,
			      dpp::module_handle_dict_type &/*moduleDict*/)
{
  // myConfig.tree_dump();

  if (myConfig.has_key("output_filename")) {
    output_filename = myConfig.fetch_string("output_filename");
    DT_LOG_DEBUG(get_logging_priority(), "setting configuration output_filename = " << output_filename);}

  if (myConfig.has_key("min_e")) {
    min_e = myConfig.fetch_integer("min_e");
    DT_LOG_DEBUG(get_logging_priority(), "setting configuration min_e = " << min_e);}
  else min_e = 0;
  
  if (myConfig.has_key("min_g")) {
    min_g = myConfig.fetch_integer("min_g");
    DT_LOG_DEBUG(get_logging_priority(), "setting configuration min_g = " << min_g);}
  else min_g = 0;
  
  ///////////////////////////////////
  
  total_entries = 0; 
  selected_entries = 0;
  
  // to allow branching vectors in tree
  gROOT->ProcessLine("#include<vector>");
  
  output_tree = new TTree ("output", "");
  output_tree->SetDirectory(0);

  // output_tree->Branch("run", &run_number);
  // output_tree->Branch("event", &event_number);

  output_tree->Branch("event",             &total_entries);

  output_tree->Branch("e", &e);
  output_tree->Branch("g", &g);

  output_tree->Branch("flag_e",                  &flag_e);
  // output_tree->Branch("ncalo_e",                &ncalo_e);
  output_tree->Branch("geomid_e",              &geomid_e);
  output_tree->Branch("energy_e",              &energy_e);
  output_tree->Branch("sigma_energy_e",  &sigma_energy_e);
  output_tree->Branch("time_e",                  &time_e);
  output_tree->Branch("sigma_time_e",      &sigma_time_e);
  // output_tree->Branch("time_th_int_e",    &time_th_int_e);
  // output_tree->Branch("time_th_ext_e",    &time_th_ext_e);
  output_tree->Branch("track_length_e",  &track_length_e);
  output_tree->Branch("x_foil_e",              &x_foil_e);
  output_tree->Branch("y_foil_e",              &y_foil_e);
  output_tree->Branch("z_foil_e",              &z_foil_e);
  output_tree->Branch("x_calo_e",              &x_calo_e); 
  output_tree->Branch("y_calo_e",              &y_calo_e);
  output_tree->Branch("z_calo_e",              &z_calo_e);
  output_tree->Branch("energy_e_true",         &energy_e_true);
  output_tree->Branch("time_e_true",             &time_e_true);
  // output_tree->Branch("x_foil_e_true",         &x_foil_e_true);
  // output_tree->Branch("y_foil_e_true",         &y_foil_e_true);
  // output_tree->Branch("z_foil_e_true",         &z_foil_e_true);
  output_tree->Branch("x_calo_e_true",         &x_calo_e_true);
  output_tree->Branch("y_calo_e_true",         &y_calo_e_true);
  output_tree->Branch("z_calo_e_true",         &z_calo_e_true);

  output_tree->Branch("flag_g",                  &flag_g);
  // output_tree->Branch("ncalo_g",                &ncalo_g);
  output_tree->Branch("geomid_g",              &geomid_g);
  output_tree->Branch("energy_g",              &energy_g);
  output_tree->Branch("sigma_energy_g",  &sigma_energy_g);
  output_tree->Branch("time_g",                  &time_g);
  output_tree->Branch("sigma_time_g",      &sigma_time_g);
  output_tree->Branch("energy_g_true",         &energy_g_true);
  output_tree->Branch("time_g_true",             &time_g_true);

  output_tree->Branch("e0", &e0);
  output_tree->Branch("g0", &g0);
  output_tree->Branch("energy_e0_true",       &energy_e0_true);
  output_tree->Branch("energy_g0_true",       &energy_g0_true);
  output_tree->Branch("x0_true",                     &x0_true);
  output_tree->Branch("y0_true",                     &y0_true);
  output_tree->Branch("z0_true",                     &z0_true);

  output_tree->Branch("pint_e1e2",   &pint_e1e2);
  output_tree->Branch("pext_e1e2",   &pext_e1e2);

  output_tree->Branch("pint_e1g1",   &pint_e1g1);
  output_tree->Branch("pext_e1g1",   &pext_e1g1);

  output_tree->Branch("pint_e2g1",   &pint_e2g1);
  output_tree->Branch("pext_e2g1",   &pext_e2g1);

  this->_set_initialized(true);
}

dpp::chain_module::process_status FLMANU_NeMg::process(datatools::things &event)
{
  if ((++total_entries % 1000000) == 0)
    DT_LOG_NOTICE(get_logging_priority(), Form("[%09lld] events", total_entries))
      
  e = 0;
  flag_e.clear();
  // ncalo_e.clear();
  geomid_e.clear();
  energy_e.clear();
  sigma_energy_e.clear();
  time_e.clear();
  sigma_time_e.clear();
  track_length_e.clear();
  time_th_int_e.clear();
  time_th_ext_e.clear();
  x_foil_e.clear();
  y_foil_e.clear();
  z_foil_e.clear();
  x_calo_e.clear();
  y_calo_e.clear();
  z_calo_e.clear();
  energy_e_true.clear();
  time_e_true.clear();
  // x_foil_e_true.clear();
  // y_foil_e_true.clear();
  // z_foil_e_true.clear();
  x_calo_e_true.clear();
  y_calo_e_true.clear();
  z_calo_e_true.clear();

  energy_e0_true.clear();
  energy_g0_true.clear();

  const snemo::datamodel::calibrated_data & CD = event.get<snemo::datamodel::calibrated_data>("CD");
  if (CD.calibrated_calorimeter_hits().size() == 0) return PROCESS_STOP;

  const snemo::datamodel::particle_track_data &PTD = event.get<snemo::datamodel::particle_track_data>("PTD");
  // if (PTD.get_number_of_particles() == 0) return PROCESS_STOP;

  mctools::simulated_data & SD = event.grab<mctools::simulated_data>("SD");

  // const genbb::primary_event &primary_event = SD.get_primary_event();
  const mctools::simulated_data::primary_event_type &primary_event = SD.get_primary_event();

  // if (!primary_event.has_vertex()) {
  //   std::cout << "*** primary event does not have vertex" << std::endl;
  //   x0_true = y0_true = z0_true = 0;}
  // else {
  //   const geomtools::vector_3d &primary_vertex = primary_event.get_vertex();
  //   x0_true = primary_vertex.x();
  //   y0_true = primary_vertex.y();
  //   z0_true = primary_vertex.z();}

  e0 = g0 = 0;

  const genbb::primary_event::particles_col_type &primary_particles = primary_event.get_particles();

  std::list<genbb::primary_particle>::const_iterator pp_it ;

  for (pp_it = primary_particles.begin(); pp_it != primary_particles.end(); ++pp_it) {
    // const std::string label = pp_it->get_particle_label();
    // const double energy = pp_it->get_kinetic_energy();
    // printf("[%06d] primary particule : '%s' E = %f MeV\n", total_entries-1, label.data(), energy);

    switch (pp_it->get_type()) {
    case genbb::primary_particle::GAMMA :
      energy_g0_true.push_back(pp_it->get_kinetic_energy());
      ++g0; break;
      
    case genbb::primary_particle::ELECTRON :
      energy_e0_true.push_back(pp_it->get_kinetic_energy());
      ++e0; break;

    case genbb::primary_particle::POSITRON :
      // energy_e0_true.push_back(pp_it->get_kinetic_energy() * -1);
      break;

    default:
      printf("[%06lld] unexpected primary particule '%s'\n", total_entries-1, pp_it->get_particle_label().data());
      break;}
  }

  x0_true = SD.get_vertex().x();
  y0_true = SD.get_vertex().y();
  z0_true = SD.get_vertex().z();
  
  
  
  
  std::vector<std::string> hit_categories;
  SD.get_step_hits_categories(hit_categories, mctools::simulated_data::HIT_CATEGORY_TYPE_PUBLIC);
  
  for (int part=0; part<PTD.get_number_of_particles(); ++part) {

    const snemo::datamodel::particle_track &part_track = PTD.get_particle(part);

    const snemo::datamodel::calibrated_calorimeter_hit::collection_type &part_hits = part_track.get_associated_calorimeter_hits();
    const snemo::datamodel::particle_track::vertex_collection_type &part_vertices = part_track.get_vertices();

    // require 1 associated calorimeter hit per track
    if (part_hits.size() != 1) return PROCESS_STOP;
    
    // require 2 vertices per track
    if (part_vertices.size() != 2) return PROCESS_STOP;
    
    const geomtools::blur_spot part_vertex0 = part_vertices[0].get();
    const geomtools::blur_spot part_vertex1 = part_vertices[1].get();

    bool part_vertex0_on_foil = snemo::datamodel::particle_track::vertex_is_on_source_foil(part_vertex0);

    bool part_vertex0_on_calo = snemo::datamodel::particle_track::vertex_is_on_main_calorimeter(part_vertex0);
    part_vertex0_on_calo |= snemo::datamodel::particle_track::vertex_is_on_x_calorimeter(part_vertex0);
    part_vertex0_on_calo |= snemo::datamodel::particle_track::vertex_is_on_gamma_veto(part_vertex0);

    bool part_vertex1_on_foil = snemo::datamodel::particle_track::vertex_is_on_source_foil(part_vertex1);

    bool part_vertex1_on_calo = snemo::datamodel::particle_track::vertex_is_on_main_calorimeter(part_vertex1);
    part_vertex1_on_calo |= snemo::datamodel::particle_track::vertex_is_on_x_calorimeter(part_vertex1);
    part_vertex1_on_calo |= snemo::datamodel::particle_track::vertex_is_on_gamma_veto(part_vertex1);
      
    // require 1 vertex on source foil and 1 vertex on (any) calo
    if ( ! ((part_vertex0_on_foil && part_vertex1_on_calo) || (part_vertex0_on_calo && part_vertex1_on_foil)) )
      return PROCESS_STOP;
    
    const geomtools::geom_id &part_hit_geom_id = part_hits[0].get().get_geom_id();

    unsigned int part_vertices_foil_id = (snemo::datamodel::particle_track::vertex_is_on_source_foil(part_vertices[0].get())) ? 0 : 1;
    unsigned int part_vertices_calo_id = (part_vertices_foil_id == 0) ? 1 : 0;

    const geomtools::blur_spot &foil_vertex = part_vertices[part_vertices_foil_id].get();
    const geomtools::blur_spot &calo_vertex = part_vertices[part_vertices_calo_id].get();


    // calorimeter hit flagging

    unsigned short _flag_e_ = 0;
    unsigned short _geomid_e_ = -1;

    switch (part_hit_geom_id.get_type()) {
     
    case 1302: // main wall
      if (part_hit_geom_id.get(1) == 0)      _flag_e_ |= mwall_it_flag;
      else if (part_hit_geom_id.get(1) == 1) _flag_e_ |= mwall_fr_flag;
      else std::cout << "*** unexpected main wall OM type = " << part_hit_geom_id.get(0) << std::endl;
      // _geomid_e_ = 10000*part_hit_geom_id.get(1) + 100*part_hit_geom_id.get(2) + part_hit_geom_id.get(3);
      _geomid_e_ = geomid_to_pmt(part_hit_geom_id);
      break;

    case 1232: // xwall
      if (part_hit_geom_id.get(1) == 0)      _flag_e_ |= xwall_it_flag;
      else if (part_hit_geom_id.get(1) == 1) _flag_e_ |= xwall_fr_flag;
      else std::cout << "*** unexpected X-wall OM type = " << part_hit_geom_id.get(0) << std::endl; 
      // _geomid_e_ = 10000*(2+part_hit_geom_id.get(1)) + 1000*part_hit_geom_id.get(2) + 100*part_hit_geom_id.get(3) + part_hit_geom_id.get(4);
      _geomid_e_ = geomid_to_pmt(part_hit_geom_id);
      break;

    case 1252: // vwall
      if (part_hit_geom_id.get(1) == 0)      _flag_e_ |= vwall_it_flag;
      else if (part_hit_geom_id.get(1) == 1) _flag_e_ |= vwall_fr_flag;
      else std::cout << "*** unexpected V-wall OM type = " << part_hit_geom_id.get(0) << std::endl;
      // _geomid_e_ = 10000*(4+part_hit_geom_id.get(1)) + 100*part_hit_geom_id.get(2) + part_hit_geom_id.get(3);
      _geomid_e_ = geomid_to_pmt(part_hit_geom_id);
      break;

    default:
      std::cout << "*** unexpected geom_id = " << part_hits[0].get().get_geom_id().get_type()<< std::endl;
      break;}
 
    float       _energy_e_ = part_hits[0].get().get_energy();
    float _sigma_energy_e_ = part_hits[0].get().get_sigma_energy() * 2.35482;

    float         _time_e_ = part_hits[0].get().get_time();
    float   _sigma_time_e_ = part_hits[0].get().get_sigma_time();

    float _x_foil_e_ = foil_vertex.get_position()[0];
    float _y_foil_e_ = foil_vertex.get_position()[1];
    float _z_foil_e_ = foil_vertex.get_position()[2];

    float _x_calo_e_ = calo_vertex.get_position()[0];
    float _y_calo_e_ = calo_vertex.get_position()[1];
    float _z_calo_e_ = calo_vertex.get_position()[2];
 
    float _track_length_e_ = 0;
    _track_length_e_ += std::pow(_x_foil_e_-_x_calo_e_, 2);
    _track_length_e_ += std::pow(_y_foil_e_-_y_calo_e_, 2);
    _track_length_e_ += std::pow(_z_foil_e_-_z_calo_e_, 2);
    _track_length_e_ = std::sqrt(_track_length_e_);


    float        _energy_e_true_ = 0;
    float          _time_e_true_ = 1E9;
    
    // float        _x_foil_e_true_ = SD.get_vertex().x();
    // float        _y_foil_e_true_ = SD.get_vertex().y();
    // float        _z_foil_e_true_ = SD.get_vertex().z();
    
    float        _x_calo_e_true_ = 0;
    float        _y_calo_e_true_ = 0;
    float        _z_calo_e_true_ = 0;

    std::vector<std::string>::const_iterator hit_category = hit_categories.begin();

    for ( ; hit_category != hit_categories.end(); ++hit_category) {
      
      const std::string &a_hit_category_name = *hit_category;
      // printf("[%06lld] particule %d has hit_category = \"%s\"\n", total_entries-1, part, a_hit_category_name.data());

      // keep only "calo" "xcalo" and "gveto"
      if (a_hit_category_name == "gg") continue;

      if (!SD.has_step_hits(a_hit_category_name)) {
	printf("[%06lld] SD does not have '%s' step\n", total_entries-1, a_hit_category_name.data());
      	continue;}
      
      mctools::simulated_data::hit_handle_collection_type &hit_collection = SD.grab_step_hits(a_hit_category_name);
    
      mctools::simulated_data::hit_handle_collection_type::iterator a_hit;

      for (a_hit = hit_collection.begin(); a_hit != hit_collection.end(); ++a_hit) {
	mctools::base_step_hit &a_step = a_hit->grab();
	
	if (a_step.get_geom_id() != part_hit_geom_id) continue;

	_energy_e_true_ += a_step.get_energy_deposit()/CLHEP::MeV;
	
	if (a_step.get_time_start()/CLHEP::ns < _time_e_true_)
	  _time_e_true_ = a_step.get_time_start()/CLHEP::ns;

	// if (_x_calo_e_true_ != 0)
	//   printf("[%06lld] multi true vertex electron\n", total_entries-1);

	// mean vertex (ponderated in energy)
	_x_calo_e_true_ += a_step.get_position_start()[0] * a_step.get_energy_deposit()/CLHEP::MeV;
	_y_calo_e_true_ += a_step.get_position_start()[1] * a_step.get_energy_deposit()/CLHEP::MeV;
	_z_calo_e_true_ += a_step.get_position_start()[2] * a_step.get_energy_deposit()/CLHEP::MeV;

      } // for a_hit

    } // for hit_category

    if (_energy_e_true_ != 0) {
      _x_calo_e_true_ /= _energy_e_true_;
      _y_calo_e_true_ /= _energy_e_true_;
      _z_calo_e_true_ /= _energy_e_true_;}
    else
      printf("[%06lld] _energy_e_true_ = 0 for particule %d !!!\n", total_entries-1, part);
      


    ++e;

    flag_e.push_back(_flag_e_);
    // ncalo_e.push_back(_ncalo_e_);

    geomid_e.push_back(_geomid_e_);

    energy_e.push_back(_energy_e_);
    sigma_energy_e.push_back(_sigma_energy_e_);

    time_e.push_back(_time_e_);
    sigma_time_e.push_back(_sigma_time_e_);

    x_foil_e.push_back(_x_foil_e_);
    y_foil_e.push_back(_y_foil_e_);
    z_foil_e.push_back(_z_foil_e_);

    x_calo_e.push_back(_x_calo_e_);
    y_calo_e.push_back(_y_calo_e_);
    z_calo_e.push_back(_z_calo_e_);

    track_length_e.push_back(_track_length_e_);


    energy_e_true.push_back(_energy_e_true_);
    time_e_true.push_back(_time_e_true_);

    // x_foil_e_true.push_back(_x_foil_e_true_);
    // y_foil_e_true.push_back(_y_foil_e_true_);
    // z_foil_e_true.push_back(_z_foil_e_true_);

    x_calo_e_true.push_back(_x_calo_e_true_);
    y_calo_e_true.push_back(_y_calo_e_true_);
    z_calo_e_true.push_back(_z_calo_e_true_);
  }

  if (e < min_e) return PROCESS_STOP;

  ///////////////////////////////////

  g = 0;
  flag_g.clear();
  // ncalo_g.clear();
  geomid_g.clear();
  energy_g.clear();
  sigma_energy_g.clear();
  time_g.clear();
  sigma_time_g.clear();

  energy_g_true.clear();
  time_g_true.clear();

  const snemo::datamodel::calibrated_calorimeter_hit::collection_type &phits = PTD.get_non_associated_calorimeters();
  // printf("[%06lld] %d unassociated hits\n", total_entries-1, phits.size());

  for (int gamma=0; gamma<phits.size(); ++gamma) {

    const geomtools::geom_id &phit_geom_id = phits[gamma].get().get_geom_id();

    unsigned short _flag_g_ = 0;
    unsigned short _geomid_g_ = -1;
    
    switch (phit_geom_id.get_type()) {
      
    case 1302: // main wall
      if (phit_geom_id.get(1) == 0)      _flag_g_ |= mwall_it_flag;
      else if (phit_geom_id.get(1) == 1) _flag_g_ |= mwall_fr_flag;
      else std::cout << "*** unexpected main wall OM type = " << phit_geom_id.get(0) << std::endl;
      // _geomid_g_ = 10000*phit_geom_id.get(1) + 100*phit_geom_id.get(2) + phit_geom_id.get(3);
      _geomid_g_ = geomid_to_pmt(phit_geom_id);
      break;
      
    case 1232: // xwall
      if (phit_geom_id.get(1) == 0)      _flag_g_ |= xwall_it_flag;
      else if (phit_geom_id.get(1) == 1) _flag_g_ |= xwall_fr_flag;
      else std::cout << "*** unexpected X-wall OM type = " << phit_geom_id.get(0) << std::endl;
      // _geomid_g_ = 10000*(2+phit_geom_id.get(1)) + 1000*phit_geom_id.get(2) + 100*phit_geom_id.get(3) + phit_geom_id.get(4);
      _geomid_g_ = geomid_to_pmt(phit_geom_id);
      break;

    case 1252: // vwall
      if (phit_geom_id.get(1) == 0)      _flag_g_ |= vwall_it_flag;
      else if (phit_geom_id.get(1) == 1) _flag_g_ |= vwall_fr_flag;
      else std::cout << "*** unexpected V-wall OM type = " << phit_geom_id.get(0) << std::endl;
      // _geomid_g_ = 10000*(4+phit_geom_id.get(1)) + 100*phit_geom_id.get(2) + phit_geom_id.get(3);
      _geomid_g_ = geomid_to_pmt(phit_geom_id);
      break;

    default:
      std::cout << "*** unexpected geom_id = " << phits[gamma].get().get_geom_id().get_type()<< std::endl;
      break;}
 
    float        _energy_g_true_ = 0;
    float          _time_g_true_ = 1E9;

    std::vector<std::string>::const_iterator hit_category = hit_categories.begin();

    for ( ; hit_category != hit_categories.end(); ++hit_category) {

      const std::string &a_hit_category_name = *hit_category;

      // keep only "calo" "xcalo" and "gveto"
      if (a_hit_category_name == "gg") continue;

      if (!SD.has_step_hits(a_hit_category_name)) {
	printf("[%06lld] SD does not have '%s' step\n", total_entries-1, a_hit_category_name.data());
      	continue;}

      mctools::simulated_data::hit_handle_collection_type &hit_collection = SD.grab_step_hits(a_hit_category_name);
    
      mctools::simulated_data::hit_handle_collection_type::iterator a_hit;

      for (a_hit = hit_collection.begin(); a_hit != hit_collection.end(); ++a_hit) {
	mctools::base_step_hit &a_step = a_hit->grab();
	
	if (a_step.get_geom_id() != phit_geom_id) continue;

	_energy_g_true_ += a_step.get_energy_deposit()/CLHEP::MeV;
	
	if (a_step.get_time_start()/CLHEP::ns < _time_g_true_)
	  _time_g_true_ = a_step.get_time_start()/CLHEP::ns;

      } // for a_hit
      
    } // for hit_category

    ++g;

    flag_g.push_back(_flag_g_);
    geomid_g.push_back(_geomid_g_);

    energy_g.push_back(phits[gamma].get().get_energy());
    sigma_energy_g.push_back(phits[gamma].get().get_sigma_energy() * 2.35482);

    time_g.push_back(phits[gamma].get().get_time());
    sigma_time_g.push_back(phits[gamma].get().get_sigma_time());

    energy_g_true.push_back(_energy_g_true_);
    time_g_true.push_back(_time_g_true_);

  }

  if (e < min_g) return PROCESS_STOP;

  ///////////////////////////////////

  pint_e1e2 = pext_e1e2 = -1;
  pint_e1g1 = pext_e1g1 = -1;
  pint_e2g1 = pext_e2g1 = -1;

  if (e == 2)
    {
      // internal probability
      
      float sigma2_time_th_int_e[2];
      
      for (int el=0; el<2; ++el) {
  	time_th_int_e.push_back(tof_tools::tof_th(energy_e[el], track_length_e[el], tof_tools::kMe));
  	sigma2_time_th_int_e[el] = tof_tools::sigma2_tof_th (time_th_int_e[el], energy_e[el], sigma_energy_e[el], tof_tools::kMe);
      }
      
      float dt_theo = time_th_int_e[1] - time_th_int_e[0];
      float dt_meas = time_e[1] - time_e[0];

      float sigma_int =
  	std::pow(sigma_time_e[0], 2)   // err meas t1
  	+ std::pow(sigma_time_e[1], 2) // err meas t2
  	+ (0.1*CLHEP::ns)*(0.1*CLHEP::ns); 
      // + sigma2_time_th_int_e1       // err theo t1
      // + sigma2_time_th_int_e2;      // err theo t2
      
      float chi2_int = std::pow(dt_meas-dt_theo, 2) / sigma_int;
        pint_e1e2 = gsl_cdf_chisq_Q (chi2_int, 1); // NDF = 1


  	// external probability

  	float sigma2_time_th_ext_e[2];

  	if (time_e[0] < time_e[1]) { //  |t1|-->--|-->--|t2| 
	  time_th_ext_e.push_back(0);
  	  time_th_ext_e.push_back(tof_tools::tof_th(energy_e[1], track_length_e[0]+track_length_e[1], tof_tools::kMe));}
  	  
  	else { //  |t1|--<--|--<--|t2| 
  	  time_th_ext_e.push_back(tof_tools::tof_th(energy_e[0], track_length_e[0]+track_length_e[1], tof_tools::kMe));
  	  time_th_ext_e.push_back(0);}
	
  	sigma2_time_th_ext_e[0] = tof_tools::sigma2_tof_th (time_th_ext_e[0], energy_e[0], sigma_energy_e[0], tof_tools::kMe);
  	sigma2_time_th_ext_e[1] = tof_tools::sigma2_tof_th (time_th_ext_e[1], energy_e[1], sigma_energy_e[1], tof_tools::kMe);

  	dt_theo = time_th_ext_e[1] - time_th_ext_e[0];
  	dt_meas = time_e[1] - time_e[0];
	
  	float sigma_ext =
  	std::pow(sigma_time_e[0], 2)   // err meas t1
  	+ std::pow(sigma_time_e[1], 2) // err meas t2
  	+ (0.1*CLHEP::ns)*(0.1*CLHEP::ns); 
  	// + sigma2_time_th_ext_e[0]       // err theo t1
  	// + sigma2_time_th_ext_e[1];      // err theo t2
  
  	float chi2_ext = std::pow(dt_meas-dt_theo, 2) / sigma_ext;
  	pext_e1e2 = gsl_cdf_chisq_Q (chi2_ext, 1); // NDF = 1
    }

  // if (g == 1) {
  //   if (e == 1) {

  //     time_th_int_e.push_back(tof_tools::tof_th(energy_e[0], track_length_e[0], tof_tools::kMe));

  //     geomtools::vector_3d vertex_e (x_foil_e, y_foil_e, z_foil_e);

  //     auto the_vertices = pte_.get_vertices();
  //     geomtools::vector_3d electron_foil_vertex;
  //     geomtools::invalidate(electron_foil_vertex);
  //     for (auto& ivtx : the_vertices) {
  //       auto a_vertex = ivtx.get();
  //       if (! snemo::datamodel::particle_track::vertex_is_on_source_foil(a_vertex))
  //         continue;

  //       electron_foil_vertex = a_vertex.get_position();
  //       break;

  //     float time_th_int_g = tof_tools::tof_th(energy_e[0], track_length_e[0], tof_tools::kMe))
  //     float dt_theo = time_th_int_e[1] - time_th_int_e[0];
  //     float dt_meas = time_e[1] - time_e[0];

  //     float sigma_int =
  // 	std::pow(sigma_time_e[0], 2)   // err meas t1
  // 	+ std::pow(sigma_time_e[1], 2) // err meas t2
  // 	+ (0.1*CLHEP::ns)*(0.1*CLHEP::ns); 
  //     // + sigma2_time_th_int_e1       // err theo t1
  //     // + sigma2_time_th_int_e2;      // err theo t2
      
  //     float chi2_int = std::pow(dt_meas-dt_theo, 2) / sigma_int;
  //       pint_e1e2 = gsl_cdf_chisq_Q (chi2_int, 1); // NDF = 1
      
  //   }
  // }

  ///////////////////////////////////

  if ((e != 0) || (g != 0)) {
    --total_entries;
    output_tree->Fill();
    ++total_entries;}

  ++selected_entries;

  return PROCESS_SUCCESS;
}

void FLMANU_NeMg::finalize()
{
  std::cout << "FLMANU_NeMg::finalize() : " << selected_entries << "/" << total_entries << " selected" << std::endl;

  TFile *output_file = new TFile(output_filename.data(), "RECREATE");
  output_file->cd(); output_tree->Write(); output_file->Close();

  this->_set_initialized(false);
}
