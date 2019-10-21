cmake_minimum_required(VERSION 3.9)

project(ZeroHitCleaner)

find_package(Falaise REQUIRED)

add_library(ZeroHitCleaner SHARED ZeroHitCleaner.cc)
target_link_libraries(ZeroHitCleaner Falaise::FalaiseModule)
configure_file("ZeroHitCleaner.conf.in" "ZeroHitCleaner.conf" @ONLY)
