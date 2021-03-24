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

install(
    TARGETS partition
)
