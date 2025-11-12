add_executable(LogPotExec source/LogPot.cpp source/Dehnen_DF.cpp)

target_include_directories(LogPotExec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../VecUtils/include ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(LogPotExec PRIVATE VecUtilsLib)


