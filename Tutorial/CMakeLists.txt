cmake_minimum_required(VERSION 3.9)

project(SkeletonModules)

find_package(Falaise REQUIRED)

add_library(SkeletonModules SHARED
  SkeletonModule.cc
  SkeletonModule_SD.cc)

target_link_libraries(SkeletonModules Falaise::FalaiseModule)
