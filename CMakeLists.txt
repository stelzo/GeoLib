cmake_minimum_required(VERSION 3.12)
project(geolib VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_TESTING "ENABLE_TESTING" OFF)

add_compile_options(-Wall -Wextra -pedantic -Werror)

# check for Mac
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
else()
    set(MAXOSX FALSE)
endif()


set(LIBRARY_SRC
    src/geo_base.cc
    src/geo_2d.cc
)

add_library(geo STATIC
    ${LIBRARY_SRC}
)

set_target_properties(geo PROPERTIES
    VERSION ${PROJECT_VERSION}
    PUBLIC_HEADER include/geo.h
)

# make install
# FRAMEWORK for Mac (dynamic and static together)
include(GNUInstallDirs)

# make uninstall
add_custom_target(uninstall
    xargs rm < install_manifest.txt
)

install(TARGETS geo 
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# files to uninstall - edit for each public header in install_manifest.in file manually
#configure_file(install_manifest.in install_manifest.txt @ONLY)

# pkg-config for easy import of the library for third parties
#configure_file(geo.pc.in geo.pc @ONLY)
#install(FILES ${CMAKE_BINARY_DIR}/geo.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

# only for shared librarys
#set_target_properties(geo PROPERTIES SOVERSION 1)

target_include_directories(geo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# only build testing when needed
if (${ENABLE_TESTING})
    add_subdirectory(tests)
    add_custom_target(test
        tests/geotest
    )
endif()
