target_sources(dns PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}/dnsanalysis.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/fieldforce.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/integrate.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/integrates.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/nonlinear.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/drive.cpp"
)

target_sources(rpo PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}/dnsanalysis.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/fieldforce.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/nonlinear.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/rpo_integrate.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/rpo_bootstrap.cpp"
)
