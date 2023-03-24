cmake_minimum_required(VERSION 3.9)

project(FLMANU)

find_package(Falaise REQUIRED)

add_library(FLMANU SHARED
            FLMANU_CD.cc
            FLMANU_Killer.cc
            FLMANU_GENBB.cc
	    FLMANU_NeMg.cc
	    FLMANU_Test.cc)
		  
target_link_libraries(FLMANU Falaise::FalaiseModule ${ROOT_LIBRARIES})

configure_file("conf/FLMANU_CD.conf.in"    "FLMANU_CD.conf" @ONLY)
configure_file("conf/FLMANU_Killer.conf.in"  "FLMANU_Killer.conf" @ONLY)
configure_file("conf/FLMANU_GENBB.conf.in" "FLMANU_GENBB.conf" @ONLY)
configure_file("conf/FLMANU_NeMg.conf.in"  "FLMANU_NeMg.conf" @ONLY)
configure_file("conf/FLMANU_Test.conf.in"  "FLMANU_Test.conf" @ONLY)
