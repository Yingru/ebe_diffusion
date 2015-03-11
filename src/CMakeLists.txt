cmake_minimum_required (VERSION 2.6)
project (diffusion)

# load Fortran compiler
enable_language (Fortran)

set (CMAKE_Fortran_FLAGS "-ffixed-line-length-none")

add_executable(${PROJECT_NAME}
  df_main.f  
  df_rnd.f  
  JetOutput-0.7.f  
  oscarnxtev.f  
  parameter_inp.f  
  rad_fnc.f  
  rad_inp.f  
  Ch_readHydro.f
  Ch_output.f  
  )

# install compiled binary
install (TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_NAME})
