add_library(partition SHARED)

target_compile_features(partition PUBLIC cxx_std_14)
set_target_properties(partition PROPERTIES
    CXX_EXTENSIONS OFF
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(partition PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_sources(partition PRIVATE
    Chunker.cc
    ChunkIndex.cc
    ChunkReducer.cc
    CmdLineUtils.cc
    ConfigStore.cc
    Csv.cc
    FileUtils.cc
    Geometry.cc
    HtmIndex.cc
    InputLines.cc
)

target_link_libraries(partition PUBLIC
    boost_program_options
    boost_thread
)


add_executable(sph-partition sph-partition.cc)
target_sources(partition PRIVATE
ObjectIndex.cc
)
target_link_libraries(sph-partition PRIVATE
partition
boost_filesystem
)


add_executable(sph-partition-matches sph-partition-matches.cc)
target_sources(sph-partition-matches PRIVATE)
target_link_libraries(sph-partition-matches PRIVATE
partition
boost_filesystem
)


install(TARGETS
    partition
    sph-partition
    sph-partition-matches
)
