cmake_minimum_required (VERSION 3.5.1)

### To use gcc/g++ on a Macintosh, you must set the Compilers
### here, not inside the project
##if(APPLE)
##       set(CMAKE_C_COMPILER "/usr/local/bin/gcc-7")
##       set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-7")
##endif()
### TODO: for now, we use CLang for Mac
###
### In order to create PALISADE's static libraries you should enable
### the BUILD_STATIC option. For that, you run "cmake .. -DBUILD_STATIC=ON".
### After having your link completed you will find static libs
### with the suffix "_static" in ./build/libs/.
### Examples: PALISADEpke_static.a, PALISADEcore_static.a, etc.
### After you run "make install" in your build directory, you can build your custom application.
### If you need your application to be linked statically, then run "cmake .. -DBUILD_STATIC=ON"

project(demo CXX)
set(CMAKE_CXX_STANDARD 11)
option( BUILD_STATIC "Set to ON to include static versions of the library" OFF)

find_package(Palisade)

set( CMAKE_CXX_FLAGS ${PALISADE_CXX_FLAGS} )

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${PALISADE_INCLUDE} )
include_directories( ${PALISADE_INCLUDE}/third-party/include )
include_directories( ${PALISADE_INCLUDE}/core )
include_directories( ${PALISADE_INCLUDE}/pke )
### add directories for other PALISADE modules as needed for your project

link_directories( ${PALISADE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
if(BUILD_STATIC)
    set( CMAKE_EXE_LINKER_FLAGS "${PALISADE_EXE_LINKER_FLAGS} -static")
    link_libraries( ${PALISADE_STATIC_LIBRARIES} )
else()
    set( CMAKE_EXE_LINKER_FLAGS ${PALISADE_EXE_LINKER_FLAGS} )
    link_libraries( ${PALISADE_SHARED_LIBRARIES} )
endif()

### ADD YOUR EXECUTABLE(s) HERE
### add_executable( EXECUTABLE-NAME SOURCES )
###
### EXAMPLE:
### add_executable( test demo-simple-example.cpp )


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
set( SOURCESBLOOM  ${CMAKE_CURRENT_SOURCE_DIR}/src/bloom.cpp ${CMAKE_CURRENT_SOURCE_DIR}/include/bloom.hpp)




add_executable( exp1stBF exp1stBF.cpp ${SOURCESBLOOM})
add_executable( exp2ndBF exp2ndBF.cpp ${SOURCESBLOOM})
add_executable( expHEPacked expHEPacked.cpp ${SOURCESBLOOM})
add_executable( expHybridPacked expHybridPacked.cpp ${SOURCESBLOOM})


add_executable( serial1stBFTemplate serial1stBFTemplate.cpp ${SOURCESBLOOM})
add_executable( serial2ndBFTemplate serial2ndBFTemplate.cpp ${SOURCESBLOOM})
add_executable( serialHETemplate serialHETemplate.cpp ${SOURCESBLOOM})
add_executable( serialHybridTemplate serialHybridTemplate.cpp ${SOURCESBLOOM})

